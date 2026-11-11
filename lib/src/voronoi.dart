import 'dart:math' as math;

import 'math.dart';
import 'mesh.dart';
import 'slicer.dart';

/// Individual convex/Voronoi shard produced by [VoronoiFractureBuilder].
final class FractureShard {
  final int index;
  final FracVec3 seedSite;
  final FractureMesh mesh;
  final FracVec3 centroid;
  final double volume;
  final FracVec3 boundsMin;
  final FracVec3 boundsMax;

  const FractureShard({
    required this.index,
    required this.seedSite,
    required this.mesh,
    required this.centroid,
    required this.volume,
    required this.boundsMin,
    required this.boundsMax,
  });
}

/// Deterministic seeded 3D Voronoi fragmentation builder using iterative bisector half-space clipping.
abstract final class VoronoiFractureBuilder {
  /// Generates [count] deterministic 3D seed sites inside `[boundsMin, boundsMax]`.
  ///
  /// When [impactPoint] is provided with [clusterBias] in `0..1`, seeds pull inward toward
  /// [impactPoint] for finer fragmentation around an impact zone.
  static List<FracVec3> generateSites({
    required FracVec3 boundsMin,
    required FracVec3 boundsMax,
    required int count,
    int seed = 42,
    FracVec3? impactPoint,
    double clusterBias = 0.0,
  }) {
    if (count <= 0) return const [];
    final rng = math.Random(seed);
    final span = boundsMax - boundsMin;
    final bias = clusterBias.clamp(0.0, 0.92);
    final sites = <FracVec3>[];

    for (var i = 0; i < count; i++) {
      // Stratified jitter avoids coincident duplicate sites.
      final rx = 0.08 + rng.nextDouble() * 0.84;
      final ry = 0.08 + rng.nextDouble() * 0.84;
      final rz = 0.08 + rng.nextDouble() * 0.84;
      var p = FracVec3(
        boundsMin.x + span.x * rx,
        boundsMin.y + span.y * ry,
        boundsMin.z + span.z * rz,
      );
      if (impactPoint != null && bias > 0.0) {
        final pull = bias * (0.35 + 0.65 * rng.nextDouble());
        p = FracVec3.lerp(p, impactPoint, pull);
      }
      sites.add(p);
    }
    return sites;
  }

  /// Shatters [mesh] into bounded 3D Voronoi cell shards.
  static List<FractureShard> fracture({
    required FractureMesh mesh,
    int cellCount = 8,
    int seed = 42,
    List<FracVec3>? sites,
    FracVec3? impactPoint,
    double clusterBias = 0.0,
    double minShardVolume = 1e-5,
  }) {
    if (mesh.isEmpty) return const [];
    final bounds = mesh.computeBounds();
    final seedSites = sites ??
        generateSites(
          boundsMin: bounds.min,
          boundsMax: bounds.max,
          count: cellCount,
          seed: seed,
          impactPoint: impactPoint,
          clusterBias: clusterBias,
        );
    if (seedSites.isEmpty) return const [];
    if (seedSites.length == 1) {
      final vol = mesh.computeSignedVolume().abs();
      return [
        FractureShard(
          index: 0,
          seedSite: seedSites.first,
          mesh: mesh,
          centroid: mesh.computeCentroid(),
          volume: vol,
          boundsMin: bounds.min,
          boundsMax: bounds.max,
        ),
      ];
    }

    final shards = <FractureShard>[];
    for (var i = 0; i < seedSites.length; i++) {
      final siteI = seedSites[i];
      var cellMesh = mesh;
      var cellBounds = bounds;

      var empty = false;
      for (var j = 0; j < seedSites.length; j++) {
        if (i == j) continue;
        final siteJ = seedSites[j];
        final delta = siteI - siteJ;
        if (delta.lengthSquared <= 1e-10) continue;
        final mid = (siteI + siteJ) * 0.5;
        final plane = FracPlane.fromPointNormal(mid, delta);

        // Quick AABB half-space rejection avoids redundant clipping passes.
        var minCorner = double.infinity;
        var maxCorner = -double.infinity;
        final bMin = cellBounds.min;
        final bMax = cellBounds.max;
        for (var c = 0; c < 8; c++) {
          final cx = (c & 1) == 0 ? bMin.x : bMax.x;
          final cy = (c & 2) == 0 ? bMin.y : bMax.y;
          final cz = (c & 4) == 0 ? bMin.z : bMax.z;
          final d = plane.signedDistance(FracVec3(cx, cy, cz));
          if (d < minCorner) minCorner = d;
          if (d > maxCorner) maxCorner = d;
        }
        if (minCorner >= -1e-5) continue;
        if (maxCorner <= 1e-5) {
          empty = true;
          break;
        }

        final res = MeshPlaneSlicer.slice(cellMesh, plane);
        if (res.positive == null || res.positive!.isEmpty) {
          empty = true;
          break;
        }
        cellMesh = res.positive!;
        cellBounds = cellMesh.computeBounds();
      }

      if (empty || cellMesh.isEmpty) continue;
      final vol = cellMesh.computeSignedVolume().abs();
      if (vol < minShardVolume) continue;

      shards.add(
        FractureShard(
          index: shards.length,
          seedSite: siteI,
          mesh: cellMesh,
          centroid: cellMesh.computeCentroid(),
          volume: vol,
          boundsMin: cellBounds.min,
          boundsMax: cellBounds.max,
        ),
      );
    }

    return shards;
  }
}
