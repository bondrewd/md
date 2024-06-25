const std = @import("std");
const F = @import("precision.zig").F;

const Coordinate = struct {
    id: u64,
    x: F,
    y: F,
    z: F,
};

pub const Crd = struct {
    coordinates: []Coordinate,

    const Self = @This();

    pub fn writeLog(self: *const Self, writer: std.fs.File.Writer) !void {
        try writer.print("[INFO] CRD\n", .{});
        try writer.print("[INFO] number of particles: {d}\n", .{self.coordinates.len});
        try writer.print("[INFO]\n", .{});
    }
};
