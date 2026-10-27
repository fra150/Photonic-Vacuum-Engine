import bpy

sc = bpy.context.scene
print("SCENARIO:", sc.name, "frames", sc.frame_start, sc.frame_end, "fps", sc.render.fps,
      "res", sc.render.resolution_x, sc.render.resolution_y)

roots = [o for o in sc.objects if "root" in o.name.lower() or "empty" in o.name.lower()]
print("ROOT-LIKE:", [(o.name, o.type) for o in roots])

objs = [o for o in sc.objects if o.type in {"MESH", "CURVE", "FONT"}]
print("OBJ COUNT:", len(objs))

for f in (1, 30, 45, 60, 75, 90, 105, 120, 150, 160, 180, 200, 220, 240):
    sc.frame_set(f)
    for o in roots:
        print("f%03d %-24s loc=(%.3f, %.3f, %.3f)" % (f, o.name, o.matrix_world.translation.x,
                                                      o.matrix_world.translation.y,
                                                      o.matrix_world.translation.z))
