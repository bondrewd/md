const std = @import("std");
const System = @import("system.zig").System;

pub const XYZFile = struct {
    file: ?std.fs.File = null,
    writer: ?std.fs.File.Writer = null,

    const Self = @This();

    pub fn init(path: ?[]const u8) Self {
        var xyz_file = Self{};
        if (path) |xyz_file_path| {
            xyz_file.file = std.fs.cwd().createFile(xyz_file_path, .{ .truncate = true }) catch {
                std.debug.print("[ERROR] Unable to open file '{s}'\n", .{xyz_file_path});
                std.process.exit(0);
            };
            xyz_file.writer = xyz_file.file.?.writer();
        }
        return xyz_file;
    }

    pub fn deinit(self: *const Self) void {
        if (self.file) |file| {
            file.close();
        }
    }

    pub fn writeFrame(self: *const Self, system: *const System) !void {
        if (self.writer) |writer| {
            try writer.print("{d}\n\n", .{system.n});
            for (system.r) |r| {
                try writer.print("Ar {d:8.3} {d:8.3} {d:8.3}\n", .{ r[0], r[1], r[2] });
            }
        }
    }
};
