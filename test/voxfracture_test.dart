import 'package:test/test.dart';
import 'package:voxfracture/voxfracture.dart';

void main() {
  group('MeshPlaneSlicer', () {
    test('slices a 2x2x2 box obliquely while conserving watertight volume and outward normals', () {
      final box = FractureMesh.box(size: const FracVec3(2.0, 2.0, 2.0));
      expect(box.computeSignedVolume(), closeTo(8.0, 1e-4));

      final plane = FracPlane.fromPointNormal(
        const FracVec3(0.0, 0.25, 0.0),
        const FracVec3(0.25, 1.0, -0.15),
      );
      final res = MeshPlaneSlicer.slice(box, plane);

      expect(res.hadCut, isTrue);
      expect(res.positive, isNotNull);
      expect(res.negative, isNotNull);

      final vPos = res.positive!.computeSignedVolume();
      final vNeg = res.negative!.computeSignedVolume();
      expect(vPos, greaterThan(1.0));
      expect(vNeg, greaterThan(2.5));
      expect(vPos + vNeg, closeTo(8.0, 1e-2));
      expect(res.positive!.computeNormalWindingAlignment(), greaterThan(0.95));
      expect(res.negative!.computeNormalWindingAlignment(), greaterThan(0.95));
    });

    test('returns hadCut == false when slice plane misses mesh completely', () {
      final box = FractureMesh.box(size: const FracVec3(1.0, 1.0, 1.0));
      final plane = FracPlane.fromPointNormal(
        const FracVec3(0.0, 5.0, 0.0),
        FracVec3.unitY,
      );
      final res = MeshPlaneSlicer.slice(box, plane);
      expect(res.hadCut, isFalse);
      expect(res.positive, isNull);
      expect(res.negative, isNotNull);
      expect(res.negative!.computeSignedVolume(), closeTo(1.0, 1e-4));
    });
  });

  group('VoronoiFractureBuilder', () {
    test('fractures box into convex shards whose volumes sum cleanly to total', () {
      final box = FractureMesh.box(size: const FracVec3(2.0, 2.0, 2.0));
      final shards = VoronoiFractureBuilder.fracture(
        mesh: box,
        cellCount: 8,
        seed: 1337,
        impactPoint: const FracVec3(0.2, 0.3, 0.0),
        clusterBias: 0.35,
      );

      expect(shards.length, greaterThanOrEqualTo(6));
      var totalVol = 0.0;
      for (final s in shards) {
        expect(s.volume, greaterThan(0.0));
        expect(s.mesh.triangleCount, greaterThanOrEqualTo(4));
        expect(s.mesh.computeNormalWindingAlignment(), greaterThan(0.9));
        totalVol += s.volume;
      }
      expect(totalVol, closeTo(8.0, 0.08));
    });
  });

  group('StructuralSupportGraph', () {
    test('detaches unreachable upper pillar pieces and over-stressed cantilevers', () {
      final graph = StructuralSupportGraph()
        ..addNode(SupportNode(id: 'base', position: const FracVec3(0, 0, 0), isAnchor: true))
        ..addNode(SupportNode(id: 'mid', position: const FracVec3(0, 1, 0), mass: 2.0))
        ..addNode(SupportNode(id: 'top', position: const FracVec3(0, 2, 0), mass: 2.0))
        ..addNode(SupportNode(id: 'overhang', position: const FracVec3(2.0, 2, 0), mass: 5.0))
        ..addEdge(SupportEdge(nodeA: 'base', nodeB: 'mid', bondStrength: 50.0))
        ..addEdge(SupportEdge(nodeA: 'mid', nodeB: 'top', bondStrength: 50.0))
        ..addEdge(SupportEdge(nodeA: 'top', nodeB: 'overhang', bondStrength: 100.0));

      final destroyMid = graph.destroyNode('mid');
      expect(destroyMid.detachedNodeIds, containsAll(['mid', 'top', 'overhang']));
      expect(destroyMid.intactCount, 1);
      expect(destroyMid.stabilityRatio, closeTo(0.25, 1e-4));
    });

    test('snaps weak horizontal cantilever bond while keeping anchor standing', () {
      final graph = StructuralSupportGraph()
        ..addNode(SupportNode(id: 'pillar', position: const FracVec3(0, 0, 0), isAnchor: true))
        ..addNode(SupportNode(id: 'arm1', position: const FracVec3(2, 0, 0), mass: 6.0))
        ..addNode(SupportNode(id: 'arm2', position: const FracVec3(4, 0, 0), mass: 8.0))
        ..addEdge(SupportEdge(nodeA: 'pillar', nodeB: 'arm1', bondStrength: 18.0))
        ..addEdge(SupportEdge(nodeA: 'arm1', nodeB: 'arm2', bondStrength: 15.0));

      final res = graph.recomputeStability();
      expect(res.brokenEdgeCount, greaterThanOrEqualTo(1));
      expect(res.detachedNodeIds, containsAll(['arm1', 'arm2']));
      expect(graph.nodes['pillar']!.intact, isTrue);
    });
  });
}
