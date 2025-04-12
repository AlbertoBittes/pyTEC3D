import vtk
import re


# js_array.SetTuple(i, js_tuple)


if __name__ == "__main__":
    filename = "3D_loop_cond_loads"
    vtk_grid = "3D_loop_cond.vtu"
    filetype = "out"
    time_steps = 5
    for i in range(1, time_steps + 1):
        element_data = parse_element_data(f"{filename}_{i}." + filetype)
        update_vtk_grid(vtk_grid, element_data, i)
