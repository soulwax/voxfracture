import 'dart:math' as math;

import 'math.dart';

/// Vertex tuple stored inside a [FractureMesh].
final class FractureVertex {
  final FracVec3 position;
  final FracVec3 normal;
  final FracVec2 uv;

  const FractureVertex({
    required this.position,
    required this.normal,
    this.uv = FracVec2.zero,
  });

  static FractureVertex lerp(FractureVertex a, FractureVertex b, double t) =>
      FractureVertex(
        position: FracVec3.lerp(a.position, b.position, t),
        normal: FracVec3.lerp(a.normal, b.normal, t).normalized,
        uv: FracVec2.lerp(a.uv, b.uv, t),
      );
}

/// Indexed triangle mesh suitable for half-space clipping, Voronoi shattering, and volume queries.
final class FractureMesh {
  final List<FractureVertex> vertices;
  final List<int> indices;

  FractureMesh({
    required this.vertices,
    required this.indices,
  });

  int get vertexCount => vertices.length;
  int get triangleCount => indices.length ~/ 3;
  bool get isEmpty => indices.length < 3 || vertices.isEmpty;

  /// Creates an axis-aligned watertight box centered at [center].
  factory FractureMesh.box({
    FracVec3 size = const FracVec3(1.0, 1.0, 1.0),
    FracVec3 center = FracVec3.zero,
  }) {
    final hx = size.x * 0.5;
    final hy = size.y * 0.5;
    final hz = size.z * 0.5;

    final verts = <FractureVertex>[];
    final idx = <int>[];

    void addQuad(
      FracVec3 p0,
      FracVec3 p1,
      FracVec3 p2,
      FracVec3 p3,
      FracVec3 n,
    ) {
      final base = verts.length;
      verts.add(FractureVertex(position: center + p0, normal: n, uv: const FracVec2(0, 0)));
      verts.add(FractureVertex(position: center + p1, normal: n, uv: const FracVec2(1, 0)));
      verts.add(FractureVertex(position: center + p2, normal: n, uv: const FracVec2(1, 1)));
      verts.add(FractureVertex(position: center + p3, normal: n, uv: const FracVec2(0, 1)));
      idx.addAll([base, base + 1, base + 2, base, base + 2, base + 3]);
    }

    // +Z
    addQuad(
      FracVec3(-hx, -hy, hz),
      FracVec3(hx, -hy, hz),
      FracVec3(hx, hy, hz),
      FracVec3(-hx, hy, hz),
      FracVec3.unitZ,
    );
    // -Z
    addQuad(
      FracVec3(hx, -hy, -hz),
      FracVec3(-hx, -hy, -hz),
      FracVec3(-hx, hy, -hz),
      FracVec3(hx, hy, -hz),
      -FracVec3.unitZ,
    );
    // +Y
    addQuad(
      FracVec3(-hx, hy, hz),
      FracVec3(hx, hy, hz),
      FracVec3(hx, hy, -hz),
      FracVec3(-hx, hy, -hz),
      FracVec3.unitY,
    );
    // -Y
    addQuad(
      FracVec3(-hx, -hy, -hz),
      FracVec3(hx, -hy, -hz),
      FracVec3(hx, -hy, hz),
      FracVec3(-hx, -hy, hz),
      -FracVec3.unitY,
    );
    // +X
    addQuad(
      FracVec3(hx, -hy, hz),
      FracVec3(hx, -hy, -hz),
      FracVec3(hx, hy, -hz),
      FracVec3(hx, hy, hz),
      FracVec3.unitX,
    );
    // -X
    addQuad(
      FracVec3(-hx, -hy, -hz),
      FracVec3(-hx, -hy, hz),
      FracVec3(-hx, hy, hz),
      FracVec3(-hx, hy, -hz),
      -FracVec3.unitX,
    );

    return FractureMesh(vertices: verts, indices: idx);
  }

  /// Creates a mesh from raw packed arrays (`positions` 3N, `normals` 3N, `uvs` 2N).
  factory FractureMesh.fromArrays({
    required List<double> positions,
    required List<double> normals,
    required List<double> uvs,
    required List<int> indices,
  }) {
    final count = positions.length ~/ 3;
    final verts = List<FractureVertex>.generate(count, (i) {
      final p3 = i * 3;
      final u2 = i * 2;
      final n = normals.length >= p3 + 3
          ? FracVec3(normals[p3], normals[p3 + 1], normals[p3 + 2]).normalized
          : FracVec3.unitY;
      final uv = uvs.length >= u2 + 2
          ? FracVec2(uvs[u2], uvs[u2 + 1])
          : FracVec2.zero;
      return FractureVertex(
        position: FracVec3(positions[p3], positions[p3 + 1], positions[p3 + 2]),
        normal: n,
        uv: uv,
      );
    });
    return FractureMesh(vertices: verts, indices: List<int>.of(indices));
  }

  /// Computes local AABB bounds `(min, max)`.
  ({FracVec3 min, FracVec3 max}) computeBounds() {
    if (vertices.isEmpty) {
      return (min: FracVec3.zero, max: FracVec3.zero);
    }
    var minX = double.infinity, minY = double.infinity, minZ = double.infinity;
    var maxX = -double.infinity, maxY = -double.infinity, maxZ = -double.infinity;
    for (final v in vertices) {
      final p = v.position;
      if (p.x < minX) minX = p.x;
      if (p.y < minY) minY = p.y;
      if (p.z < minZ) minZ = p.z;
      if (p.x > maxX) maxX = p.x;
      if (p.y > maxY) maxY = p.y;
      if (p.z > maxZ) maxZ = p.z;
    }
    return (
      min: FracVec3(minX, minY, minZ),
      max: FracVec3(maxX, maxY, maxZ),
    );
  }

  /// Computes exact signed volume via the Divergence Theorem over oriented triangles.
  double computeSignedVolume() {
    if (isEmpty) return 0.0;
    var sum = 0.0;
    for (var i = 0; i + 2 < indices.length; i += 3) {
      final a = vertices[indices[i]].position;
      final b = vertices[indices[i + 1]].position;
      final c = vertices[indices[i + 2]].position;
      sum += a.dot(b.cross(c));
    }
    return sum / 6.0;
  }

  /// Computes volume-weighted centroid (falls back to vertex centroid for zero-volume inputs).
  FracVec3 computeCentroid() {
    if (vertices.isEmpty) return FracVec3.zero;
    var volSum = 0.0;
    var weightedSum = FracVec3.zero;
    for (var i = 0; i + 2 < indices.length; i += 3) {
      final a = vertices[indices[i]].position;
      final b = vertices[indices[i + 1]].position;
      final c = vertices[indices[i + 2]].position;
      final tet6 = a.dot(b.cross(c));
      volSum += tet6;
      weightedSum += (a + b + c) * (tet6 * 0.25);
    }
    if (volSum.abs() > 1e-9) {
      return weightedSum / volSum;
    }
    var acc = FracVec3.zero;
    for (final v in vertices) {
      acc += v.position;
    }
    return acc / vertices.length.toDouble();
  }

  /// Returns a translated copy shifted by [offset] (e.g. centering a shard around its centroid).
  FractureMesh translated(FracVec3 offset) => FractureMesh(
        vertices: [
          for (final v in vertices)
            FractureVertex(
              position: v.position + offset,
              normal: v.normal,
              uv: v.uv,
            ),
        ],
        indices: List<int>.of(indices),
      );

  /// Computes area-weighted average face/vertex normal alignment sanity score `[-1..1]`.
  double computeNormalWindingAlignment() {
    var weightSum = 0.0;
    var dotSum = 0.0;
    for (var i = 0; i + 2 < indices.length; i += 3) {
      final v0 = vertices[indices[i]];
      final v1 = vertices[indices[i + 1]];
      final v2 = vertices[indices[i + 2]];
      final cross = (v1.position - v0.position).cross(v2.position - v0.position);
      final area2 = cross.length;
      if (area2 <= 1e-10) continue;
      final faceN = cross / area2;
      final avgStored = (v0.normal + v1.normal + v2.normal).normalized;
      dotSum += faceN.dot(avgStored) * area2;
      weightSum += area2;
    }
    return weightSum > 1e-12 ? dotSum / weightSum : 1.0;
  }
}

/// Utility helper returning whether two bounding boxes overlap on all axes.
bool fracBoundsOverlap(
  FracVec3 aMin,
  FracVec3 aMax,
  FracVec3 bMin,
  FracVec3 bMax, {
  double padding = 0.0,
}) =>
    aMin.x <= bMax.x + padding &&
    aMax.x >= bMin.x - padding &&
    aMin.y <= bMax.y + padding &&
    aMax.y >= bMin.y - padding &&
    aMin.z <= bMax.z + padding &&
    aMax.z >= bMin.z - padding &&
    math.min(aMax.x - aMin.x, 1.0) >= 0.0;
