import 'dart:collection';
import 'dart:math' as math;

import 'math.dart';

/// Structural element inside a [StructuralSupportGraph].
final class SupportNode {
  final String id;
  final FracVec3 position;
  final double mass;
  bool isAnchor;
  double hitPoints;
  final double maxHitPoints;
  bool intact;

  SupportNode({
    required this.id,
    required this.position,
    this.mass = 1.0,
    this.isAnchor = false,
    double maxHitPoints = 100.0,
    this.intact = true,
  })  : maxHitPoints = math.max(1.0, maxHitPoints),
        hitPoints = math.max(1.0, maxHitPoints);
}

/// Load-bearing structural connection between two [SupportNode]s.
final class SupportEdge {
  final String nodeA;
  final String nodeB;

  /// Maximum accumulated structural force/torque load supported across this joint.
  final double bondStrength;
  bool intact;

  SupportEdge({
    required this.nodeA,
    required this.nodeB,
    this.bondStrength = 24.0,
    this.intact = true,
  });

  String other(String id) => id == nodeA ? nodeB : nodeA;
}

/// Summary produced whenever damage or structural updates alter support connectivity.
final class CollapseResult {
  /// IDs of nodes that lost support in this update (destroyed + newly collapsed).
  final List<String> detachedNodeIds;

  /// Edges broken due to excessive cantilever or mass load.
  final int brokenEdgeCount;

  /// Remaining intact anchored nodes.
  final int intactCount;

  /// Fraction of original structural nodes still standing (`0..1`).
  final double stabilityRatio;

  const CollapseResult({
    required this.detachedNodeIds,
    required this.brokenEdgeCount,
    required this.intactCount,
    required this.stabilityRatio,
  });
}

/// Tracks grounded anchor reachability and cantilever load limits for destructible structures.
final class StructuralSupportGraph {
  final Map<String, SupportNode> _nodes = {};
  final List<SupportEdge> _edges = [];
  final Map<String, List<SupportEdge>> _adjacency = {};

  Map<String, SupportNode> get nodes => Map.unmodifiable(_nodes);
  List<SupportEdge> get edges => List.unmodifiable(_edges);

  int get intactCount {
    var count = 0;
    for (final n in _nodes.values) {
      if (n.intact) count++;
    }
    return count;
  }

  double get stabilityRatio =>
      _nodes.isEmpty ? 1.0 : intactCount / _nodes.length.toDouble();

  void addNode(SupportNode node) {
    _nodes[node.id] = node;
    _adjacency.putIfAbsent(node.id, () => []);
  }

  void addEdge(SupportEdge edge) {
    if (!_nodes.containsKey(edge.nodeA) || !_nodes.containsKey(edge.nodeB)) {
      throw ArgumentError('Both endpoints (${edge.nodeA}, ${edge.nodeB}) must exist in graph.');
    }
    _edges.add(edge);
    _adjacency[edge.nodeA]!.add(edge);
    _adjacency[edge.nodeB]!.add(edge);
  }

  /// Destroys [nodeId] directly and collapses any unsupported islands or overloaded cantilevers.
  CollapseResult destroyNode(String nodeId) {
    final node = _nodes[nodeId];
    if (node == null || !node.intact) {
      return recomputeStability();
    }
    node.intact = false;
    node.hitPoints = 0.0;
    return recomputeStability(newlyDestroyed: [nodeId]);
  }

  /// Applies spherical blast/impact damage centered at [impactPoint] and returns collapsed nodes.
  CollapseResult applyDamageAt(
    FracVec3 impactPoint, {
    double radius = 2.0,
    double damage = 150.0,
  }) {
    if (radius <= 0.0 || damage <= 0.0) {
      return recomputeStability();
    }
    final destroyed = <String>[];
    for (final node in _nodes.values) {
      if (!node.intact) continue;
      final dist = node.position.distanceTo(impactPoint);
      if (dist > radius) continue;
      final falloff = 1.0 - (dist / radius);
      node.hitPoints -= damage * falloff;
      if (node.hitPoints <= 0.0) {
        node.hitPoints = 0.0;
        node.intact = false;
        destroyed.add(node.id);
      }
    }
    return recomputeStability(newlyDestroyed: destroyed);
  }

  /// Propagates support from all intact anchored nodes, snapping overloaded cantilever bonds
  /// and detaching all unreachable floating islands.
  CollapseResult recomputeStability({
    List<String> newlyDestroyed = const [],
    double cantileverLeverageFactor = 1.5,
  }) {
    final detached = <String>[...newlyDestroyed];
    var brokenEdges = 0;

    while (true) {
      final parentEdge = <String, SupportEdge?>{};
      final parentId = <String, String?>{};
      final bfsOrder = <String>[];
      final queue = Queue<String>();

      for (final node in _nodes.values) {
        if (node.intact && node.isAnchor) {
          parentId[node.id] = null;
          parentEdge[node.id] = null;
          queue.add(node.id);
          bfsOrder.add(node.id);
        }
      }

      while (queue.isNotEmpty) {
        final currId = queue.removeFirst();
        for (final edge in _adjacency[currId] ?? const <SupportEdge>[]) {
          if (!edge.intact) continue;
          final nextId = edge.other(currId);
          final nextNode = _nodes[nextId];
          if (nextNode == null || !nextNode.intact) continue;
          if (parentId.containsKey(nextId)) continue;
          parentId[nextId] = currId;
          parentEdge[nextId] = edge;
          queue.add(nextId);
          bfsOrder.add(nextId);
        }
      }

      // Accumulate structural load bottom-up from leaf nodes toward anchors.
      final subtreeLoad = <String, double>{
        for (final id in bfsOrder) id: _nodes[id]!.mass,
      };
      var fracturedBondThisPass = false;
      for (var i = bfsOrder.length - 1; i >= 0; i--) {
        final id = bfsOrder[i];
        final pId = parentId[id];
        final pEdge = parentEdge[id];
        if (pId == null || pEdge == null) continue;

        final nodePos = _nodes[id]!.position;
        final parentPos = _nodes[pId]!.position;
        final dx = nodePos.x - parentPos.x;
        final dz = nodePos.z - parentPos.z;
        final horizontalSpan = math.sqrt(dx * dx + dz * dz);
        final loadThroughEdge =
            subtreeLoad[id]! * (1.0 + horizontalSpan * cantileverLeverageFactor);

        if (loadThroughEdge > pEdge.bondStrength) {
          pEdge.intact = false;
          brokenEdges++;
          fracturedBondThisPass = true;
        } else {
          subtreeLoad[pId] = (subtreeLoad[pId] ?? 0.0) + subtreeLoad[id]!;
        }
      }

      if (fracturedBondThisPass) continue;

      // Any intact node not reached from a live anchor collapses.
      for (final node in _nodes.values) {
        if (node.intact && !parentId.containsKey(node.id)) {
          node.intact = false;
          detached.add(node.id);
        }
      }
      break;
    }

    return CollapseResult(
      detachedNodeIds: detached,
      brokenEdgeCount: brokenEdges,
      intactCount: intactCount,
      stabilityRatio: stabilityRatio,
    );
  }
}
