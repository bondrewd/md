const std = @import("std");
const System = @import("system.zig").System;

pub const CSVFile = struct {
    file: ?std.fs.File = null,
    writer: ?std.fs.File.Writer = null,

    const Self = @This();

    pub fn init(path: ?[]const u8) Self {
        var csv_file = Self{};
        if (path) |csv_file_path| {
            csv_file.file = std.fs.cwd().createFile(csv_file_path, .{ .truncate = true }) catch {
                std.debug.print("[ERROR] Unable to open file '{s}'\n", .{csv_file_path});
                std.process.exit(0);
            };
            csv_file.writer = csv_file.file.?.writer();
        }
        return csv_file;
    }

    pub fn deinit(self: *const Self) void {
        if (self.file) |file| {
            file.close();
        }
    }

    pub fn writeHeader(self: *const Self) !void {
        if (self.writer) |writer| {
            try writer.print("{s},{s},{s},{s},{s},{s}\n", .{
                "STEP",
                "TEMP",
                "TOTAL",
                "POTENTIAL",
                "KINETIC",
                "LJ",
            });
        }
    }

    pub fn writeState(self: *const Self, system: *const System, step: usize) !void {
        if (self.writer) |writer| {
            try writer.print("{d},{d:.2},{d:.4},{d:.4},{d:.4},{d:.4}\n", .{
                step,
                system.measureTemperature(),
                system.energy.lj + system.energy.kinetic,
                system.energy.lj,
                system.energy.kinetic,
                system.energy.lj,
            });
        }
    }
};
