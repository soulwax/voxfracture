import 'dart:math' as math;

/// Immutable 2D vector used for UV coordinates and cap polygon triangulation.
final class FracVec2 {
  final double x;
  final double y;

  const FracVec2(this.x, this.y);

  static const FracVec2 zero = FracVec2(0.0, 0.0);

  FracVec2 operator +(FracVec2 o) => FracVec2(x + o.x, y + o.y);
  FracVec2 operator -(FracVec2 o) => FracVec2(x - o.x, y - o.y);
  FracVec2 operator *(double s) => FracVec2(x * s, y * s);

  double cross(FracVec2 o) => x * o.y - y * o.x;
  double dot(FracVec2 o) => x * o.x + y * o.y;

  static FracVec2 lerp(FracVec2 a, FracVec2 b, double t) =>
      FracVec2(a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t);

  @override
  String toString() => 'FracVec2($x, $y)';
}

/// Immutable 3D vector for [voxfracture].
final class FracVec3 {
  final double x;
  final double y;
  final double z;

  const FracVec3(this.x, this.y, this.z);

  static const FracVec3 zero = FracVec3(0.0, 0.0, 0.0);
  static const FracVec3 unitX = FracVec3(1.0, 0.0, 0.0);
  static const FracVec3 unitY = FracVec3(0.0, 1.0, 0.0);
  static const FracVec3 unitZ = FracVec3(0.0, 0.0, 1.0);

  FracVec3 operator +(FracVec3 o) => FracVec3(x + o.x, y + o.y, z + o.z);
  FracVec3 operator -(FracVec3 o) => FracVec3(x - o.x, y - o.y, z - o.z);
  FracVec3 operator -() => FracVec3(-x, -y, -z);
  FracVec3 operator *(double s) => FracVec3(x * s, y * s, z * s);
  FracVec3 operator /(double s) => FracVec3(x / s, y / s, z / s);

  double dot(FracVec3 o) => x * o.x + y * o.y + z * o.z;

  FracVec3 cross(FracVec3 o) => FracVec3(
        y * o.z - z * o.y,
        z * o.x - x * o.z,
        x * o.y - y * o.x,
      );

  double get lengthSquared => x * x + y * y + z * z;
  double get length => math.sqrt(lengthSquared);

  FracVec3 get normalized {
    final lenSq = lengthSquared;
    if (lenSq <= 1e-16) return unitY;
    final inv = 1.0 / math.sqrt(lenSq);
    return FracVec3(x * inv, y * inv, z * inv);
  }

  double distanceTo(FracVec3 o) => (this - o).length;

  static FracVec3 lerp(FracVec3 a, FracVec3 b, double t) => FracVec3(
        a.x + (b.x - a.x) * t,
        a.y + (b.y - a.y) * t,
        a.z + (b.z - a.z) * t,
      );

  @override
  String toString() => 'FracVec3($x, $y, $z)';
}

/// 3D slicing plane defined by `normal.dot(p) - distance = 0`.
final class FracPlane {
  final FracVec3 normal;
  final double distance;

  FracPlane({
    required FracVec3 normal,
    required this.distance,
  }) : normal = normal.normalized;

  factory FracPlane.fromPointNormal(FracVec3 point, FracVec3 normal) {
    final unit = normal.normalized;
    return FracPlane(normal: unit, distance: unit.dot(point));
  }

  double signedDistance(FracVec3 p) => normal.dot(p) - distance;

  FracPlane get flipped => FracPlane(normal: -normal, distance: -distance);
}
