import 'dart:math' as math;

import 'math.dart';
import 'mesh.dart';

/// Output of slicing a [FractureMesh] with a [FracPlane].
final class SliceResult {
  /// Geometry in the positive half-space (`plane.signedDistance(p) >= 0`).
  final FractureMesh? positive;

  /// Geometry in the negative half-space (`plane.signedDistance(p) <= 0`).
  final FractureMesh? negative;

  /// True when the plane actually intersected the input mesh on both sides.
  final bool hadCut;

  const SliceResult({
    this.positive,
    this.negative,
    required this.hadCut,
  });
}

final class _MeshAccumulator {
  final List<FractureVertex> vertices = [];
  final List<int> indices = [];
  bool hasOffPlane = false;

  void addTriangle(
    FractureVertex a,
    FractureVertex b,
    FractureVertex c, {
    required bool offPlane,
  }) {
    final cross = (b.position - a.position).cross(c.position - a.position);
    if (cross.lengthSquared <= 1e-14) return;
    if (offPlane) hasOffPlane = true;
    final base = vertices.length;
    vertices.add(a);
    vertices.add(b);
    vertices.add(c);
    indices.addAll([base, base + 1, base + 2]);
  }

  FractureMesh toMesh() => FractureMesh(vertices: vertices, indices: indices);
}

/// Exact CSG half-space mesh slicer with manifold boundary-loop ear-clipping caps.
abstract final class MeshPlaneSlicer {
  /// Splits [mesh] along [plane] into watertight [positive] and [negative] halves.
  static SliceResult slice(
    FractureMesh mesh,
    FracPlane plane, {
    double epsilon = 1e-5,
    bool capHoles = true,
  }) {
    if (mesh.isEmpty) {
      return const SliceResult(hadCut: false);
    }

    final pos = _MeshAccumulator();
    final neg = _MeshAccumulator();
    var anyPos = false;
    var anyNeg = false;

    for (var i = 0; i + 2 < mesh.indices.length; i += 3) {
      final v0 = mesh.vertices[mesh.indices[i]];
      final v1 = mesh.vertices[mesh.indices[i + 1]];
      final v2 = mesh.vertices[mesh.indices[i + 2]];

      final d0 = plane.signedDistance(v0.position);
      final d1 = plane.signedDistance(v1.position);
      final d2 = plane.signedDistance(v2.position);

      final s0 = _sign(d0, epsilon);
      final s1 = _sign(d1, epsilon);
      final s2 = _sign(d2, epsilon);

      if (s0 > 0 || s1 > 0 || s2 > 0) anyPos = true;
      if (s0 < 0 || s1 < 0 || s2 < 0) anyNeg = true;

      // Strictly or non-strictly positive (`all >= 0`, not all 0).
      if (s0 >= 0 && s1 >= 0 && s2 >= 0 && (s0 > 0 || s1 > 0 || s2 > 0)) {
        pos.addTriangle(v0, v1, v2, offPlane: true);
        continue;
      }

      // Strictly or non-strictly negative (`all <= 0`, not all 0).
      if (s0 <= 0 && s1 <= 0 && s2 <= 0 && (s0 < 0 || s1 < 0 || s2 < 0)) {
        neg.addTriangle(v0, v1, v2, offPlane: true);
        continue;
      }

      // Coplanar (`s0 == 0 && s1 == 0 && s2 == 0`).
      if (s0 == 0 && s1 == 0 && s2 == 0) {
        final faceNormal = (v1.position - v0.position).cross(v2.position - v0.position);
        if (faceNormal.dot(plane.normal) >= 0.0) {
          neg.addTriangle(v0, v1, v2, offPlane: false);
        } else {
          pos.addTriangle(v0, v1, v2, offPlane: false);
        }
        continue;
      }

      // Straddling triangle: clip against `signedDistance >= 0` for positive,
      // and `signedDistance <= 0` for negative using exact edge interpolation.
      final polyPos = _clipTriangleToHalfSpace([v0, v1, v2], [d0, d1, d2], keepPositive: true, epsilon: epsilon);
      _emitFan(pos, polyPos);

      final polyNeg = _clipTriangleToHalfSpace([v0, v1, v2], [d0, d1, d2], keepPositive: false, epsilon: epsilon);
      _emitFan(neg, polyNeg);
    }

    if (!anyPos) {
      return SliceResult(
        positive: null,
        negative: neg.hasOffPlane ? neg.toMesh() : mesh,
        hadCut: false,
      );
    }
    if (!anyNeg) {
      return SliceResult(
        positive: pos.hasOffPlane ? pos.toMesh() : mesh,
        negative: null,
        hadCut: false,
      );
    }

    if (capHoles) {
      _capOpenCutBoundary(pos, plane, capNormal: -plane.normal, epsilon: epsilon);
      _capOpenCutBoundary(neg, plane, capNormal: plane.normal, epsilon: epsilon);
    }

    final posMesh = pos.toMesh();
    final negMesh = neg.toMesh();
    return SliceResult(
      positive: posMesh.isEmpty ? null : posMesh,
      negative: negMesh.isEmpty ? null : negMesh,
      hadCut: true,
    );
  }

  static int _sign(double d, double eps) {
    if (d > eps) return 1;
    if (d < -eps) return -1;
    return 0;
  }

  /// Sutherland-Hodgman single half-space triangle clipper preserving winding order.
  static List<FractureVertex> _clipTriangleToHalfSpace(
    List<FractureVertex> verts,
    List<double> dists, {
    required bool keepPositive,
    required double epsilon,
  }) {
    final result = <FractureVertex>[];
    for (var i = 0; i < 3; i++) {
      final j = (i + 1) % 3;
      final a = verts[i];
      final b = verts[j];
      final da = keepPositive ? dists[i] : -dists[i];
      final db = keepPositive ? dists[j] : -dists[j];
      final inA = da >= -epsilon;
      final inB = db >= -epsilon;

      if (inA && inB) {
        result.add(b);
      } else if (inA && !inB) {
        final denom = da - db;
        final t = denom.abs() > 1e-12 ? (da / denom).clamp(0.0, 1.0) : 0.0;
        result.add(FractureVertex.lerp(a, b, t));
      } else if (!inA && inB) {
        final denom = da - db;
        final t = denom.abs() > 1e-12 ? (da / denom).clamp(0.0, 1.0) : 0.0;
        result.add(FractureVertex.lerp(a, b, t));
        result.add(b);
      }
    }
    return result;
  }

  static void _emitFan(_MeshAccumulator acc, List<FractureVertex> poly) {
    if (poly.length < 3) return;
    for (var i = 1; i + 1 < poly.length; i++) {
      acc.addTriangle(poly[0], poly[i], poly[i + 1], offPlane: true);
    }
  }

  static String _quantizeKey(FracVec3 p) {
    const scale = 10000.0; // 0.1 mm weld grid for topological boundary loop stitching
    final qx = (p.x * scale).round();
    final qy = (p.y * scale).round();
    final qz = (p.z * scale).round();
    return '$qx,$qy,$qz';
  }

  /// Finds all unpaired triangle half-edges lying ON [plane] and triangulates their closed loops
  /// using 2D ear-clipping with [capNormal] outward normal.
  static void _capOpenCutBoundary(
    _MeshAccumulator acc,
    FracPlane plane, {
    required FracVec3 capNormal,
    required double epsilon,
  }) {
    final pointTol = math.max(epsilon * 4.0, 2e-4);
    final edgeCounts = <String, int>{};
    final edgeData = <String, ({String kA, String kB, FracVec3 a, FracVec3 b})>{};
    final keyPositions = <String, FracVec3>{};

    void recordEdge(FracVec3 u, FracVec3 v) {
      if ((v - u).lengthSquared <= 1e-12) return;
      final ku = _quantizeKey(u);
      final kv = _quantizeKey(v);
      if (ku == kv) return;
      keyPositions.putIfAbsent(ku, () => u);
      keyPositions.putIfAbsent(kv, () => v);
      final fwdKey = '$ku->$kv';
      final revKey = '$kv->$ku';
      final revCount = edgeCounts[revKey] ?? 0;
      if (revCount > 0) {
        if (revCount == 1) {
          edgeCounts.remove(revKey);
          edgeData.remove(revKey);
        } else {
          edgeCounts[revKey] = revCount - 1;
        }
      } else {
        edgeCounts[fwdKey] = (edgeCounts[fwdKey] ?? 0) + 1;
        edgeData[fwdKey] = (kA: ku, kB: kv, a: u, b: v);
      }
    }

    for (var i = 0; i + 2 < acc.indices.length; i += 3) {
      final p0 = acc.vertices[acc.indices[i]].position;
      final p1 = acc.vertices[acc.indices[i + 1]].position;
      final p2 = acc.vertices[acc.indices[i + 2]].position;
      final on0 = plane.signedDistance(p0).abs() <= pointTol;
      final on1 = plane.signedDistance(p1).abs() <= pointTol;
      final on2 = plane.signedDistance(p2).abs() <= pointTol;
      if (on0 && on1) recordEdge(p0, p1);
      if (on1 && on2) recordEdge(p1, p2);
      if (on2 && on0) recordEdge(p2, p0);
    }

    if (edgeData.isEmpty) return;

    // Reverse every remaining unpaired half-edge `a -> b` into cap perimeter edge `b -> a`
    // so shared seam edges run in opposite manifold directions.
    final adjacency = <String, List<String>>{};
    for (final e in edgeData.values) {
      adjacency.putIfAbsent(e.kB, () => []).add(e.kA);
    }

    final loops = <List<FracVec3>>[];
    while (adjacency.isNotEmpty) {
      final startKey = adjacency.keys.first;
      final chainKeys = <String>[];
      var curr = startKey;
      while (true) {
        final nextList = adjacency[curr];
        if (nextList == null || nextList.isEmpty) break;
        chainKeys.add(curr);
        final next = nextList.removeAt(0);
        if (nextList.isEmpty) adjacency.remove(curr);
        curr = next;
        if (curr == startKey) break;
        if (chainKeys.length > 4096) break;
      }
      if (chainKeys.length >= 3) {
        loops.add([for (final k in chainKeys) keyPositions[k]!]);
      }
    }

    if (loops.isEmpty) return;

    // Construct orthonormal tangent basis `(tU, tV)` such that `tU x tV == capNormal`.
    final ref = capNormal.x.abs() < 0.85 ? FracVec3.unitX : FracVec3.unitY;
    final tU = ref.cross(capNormal).normalized;
    final tV = capNormal.cross(tU).normalized;

    for (final rawLoop in loops) {
      final cleaned = _removeCollinearLoopPoints(rawLoop, tU, tV);
      if (cleaned.length < 3) continue;
      _triangulateCapPolygon(acc, cleaned, capNormal, tU, tV);
    }
  }

  static List<FracVec3> _removeCollinearLoopPoints(
    List<FracVec3> loop,
    FracVec3 tU,
    FracVec3 tV,
  ) {
    var pts = List<FracVec3>.of(loop);
    var changed = true;
    while (changed && pts.length > 3) {
      changed = false;
      final next = <FracVec3>[];
      final n = pts.length;
      for (var i = 0; i < n; i++) {
        final prev = pts[(i - 1 + n) % n];
        final curr = pts[i];
        final nxt = pts[(i + 1) % n];
        if ((curr - prev).lengthSquared <= 1e-12) {
          changed = true;
          continue;
        }
        final a2 = FracVec2((curr - prev).dot(tU), (curr - prev).dot(tV));
        final b2 = FracVec2((nxt - curr).dot(tU), (nxt - curr).dot(tV));
        final denom = math.sqrt(a2.dot(a2) * b2.dot(b2));
        if (denom > 1e-12 && a2.cross(b2).abs() / denom < 1e-4 && a2.dot(b2) > 0) {
          changed = true;
          continue;
        }
        next.add(curr);
      }
      pts = next;
    }
    return pts;
  }

  static void _triangulateCapPolygon(
    _MeshAccumulator acc,
    List<FracVec3> loop,
    FracVec3 capNormal,
    FracVec3 tU,
    FracVec3 tV,
  ) {
    final pts2D = [
      for (final p in loop) FracVec2(p.dot(tU), p.dot(tV)),
    ];

    // Bounds for normalized planar cap UVs.
    var minU = double.infinity, maxU = -double.infinity;
    var minV = double.infinity, maxV = -double.infinity;
    for (final q in pts2D) {
      if (q.x < minU) minU = q.x;
      if (q.x > maxU) maxU = q.x;
      if (q.y < minV) minV = q.y;
      if (q.y > maxV) maxV = q.y;
    }
    final spanU = math.max(1e-4, maxU - minU);
    final spanV = math.max(1e-4, maxV - minV);

    // Ensure CCW ordering in `(tU, tV)` (positive signed area => normal == tU x tV == capNormal).
    var signedArea2 = 0.0;
    for (var i = 0; i < pts2D.length; i++) {
      final a = pts2D[i];
      final b = pts2D[(i + 1) % pts2D.length];
      signedArea2 += a.cross(b);
    }

    final indicesList = signedArea2 >= 0
        ? [for (var i = 0; i < loop.length; i++) i]
        : [for (var i = loop.length - 1; i >= 0; i--) i];

    FractureVertex makeCapVertex(int idx) {
      final p = loop[idx];
      final q = pts2D[idx];
      return FractureVertex(
        position: p,
        normal: capNormal,
        uv: FracVec2((q.x - minU) / spanU, (q.y - minV) / spanV),
      );
    }

    // Ear clipping over simple loop indices.
    final remaining = List<int>.of(indicesList);
    var guard = 0;
    while (remaining.length > 3 && guard < remaining.length * remaining.length + 16) {
      guard++;
      var clipped = false;
      final count = remaining.length;
      for (var i = 0; i < count; i++) {
        final iPrev = remaining[(i - 1 + count) % count];
        final iCurr = remaining[i];
        final iNext = remaining[(i + 1) % count];
        final a = pts2D[iPrev];
        final b = pts2D[iCurr];
        final c = pts2D[iNext];
        if ((b - a).cross(c - a) <= 1e-10) continue;

        var emptyInside = true;
        for (var j = 0; j < count; j++) {
          final idx = remaining[j];
          if (idx == iPrev || idx == iCurr || idx == iNext) continue;
          if (_pointInTri2D(pts2D[idx], a, b, c)) {
            emptyInside = false;
            break;
          }
        }
        if (!emptyInside) continue;

        acc.addTriangle(
          makeCapVertex(iPrev),
          makeCapVertex(iCurr),
          makeCapVertex(iNext),
          offPlane: false,
        );
        remaining.removeAt(i);
        clipped = true;
        break;
      }
      if (!clipped) break;
    }

    if (remaining.length == 3) {
      acc.addTriangle(
        makeCapVertex(remaining[0]),
        makeCapVertex(remaining[1]),
        makeCapVertex(remaining[2]),
        offPlane: false,
      );
    } else if (remaining.length > 3) {
      // Fallback convex fan if numerical degeneracy prevented ear removal.
      for (var i = 1; i + 1 < remaining.length; i++) {
        acc.addTriangle(
          makeCapVertex(remaining[0]),
          makeCapVertex(remaining[i]),
          makeCapVertex(remaining[i + 1]),
          offPlane: false,
        );
      }
    }
  }

  static bool _pointInTri2D(FracVec2 p, FracVec2 a, FracVec2 b, FracVec2 c) {
    final c1 = (b - a).cross(p - a);
    final c2 = (c - b).cross(p - b);
    final c3 = (a - c).cross(p - c);
    return c1 >= -1e-9 && c2 >= -1e-9 && c3 >= -1e-9;
  }
}
