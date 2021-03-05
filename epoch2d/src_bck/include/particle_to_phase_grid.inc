        cell_x_r = (photon_px(ii) - min_px) / dpx - 0.5_num
        cell_y_r = (photon_py(ii) - min_py) / dpy - 0.5_num
        cell_z_r = (photon_pz(ii) - min_pz) / dpz - 0.5_num

        cell_x = FLOOR(cell_x_r + 0.5_num)
        cell_y = FLOOR(cell_y_r + 0.5_num)
        cell_z = FLOOR(cell_z_r + 0.5_num)

        cell_frac_x = REAL(cell_x, num) - cell_x_r
        cell_frac_y = REAL(cell_y, num) - cell_y_r
        cell_frac_z = REAL(cell_z, num) - cell_z_r

        cell_x = cell_x + 1
        cell_y = cell_y + 1
        cell_z = cell_z + 1

        cf2 = cell_frac_x**2
        gx(-1) = 0.25_num + cf2 + cell_frac_x
        gx( 0) = 1.5_num - 2.0_num * cf2
        gx( 1) = 0.25_num + cf2 - cell_frac_x

        cf2 = cell_frac_y**2
        gy(-1) = 0.25_num + cf2 + cell_frac_y
        gy( 0) = 1.5_num - 2.0_num * cf2
        gy( 1) = 0.25_num + cf2 - cell_frac_y

        cf2 = cell_frac_z**2
        gz(-1) = 0.25_num + cf2 + cell_frac_z
        gz( 0) = 1.5_num - 2.0_num * cf2
        gz( 1) = 0.25_num + cf2 - cell_frac_z


