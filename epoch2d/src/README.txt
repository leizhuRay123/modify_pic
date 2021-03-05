update at 20201008
add z_min_boundary in injector.F90
ASSOCIATED(injector_z_min) is F
inject_in_z in species_list
working_injector in boundary.f90
create_working_injector in injector.F90
as_boundary_ in deck/string.F90
calculate_density and temperature use the same function
note return behind injector_in_z_direction

add density_profile by using populate_properties of injector and also temperature and drift


to control the injection velocity, we need to understand the injector%depth similar to the optical depth

assump dz = 1 m

some error in boundary because of the differences of x_min y_min and z_min 

appears a problem that there are many particles' weight = 0, should set density_min


#########
open boundary has some problems:
open = simple_outflow
field_clamp_zero(ex,ey,ez)
field_clamp_zero(bx,by,bz)

Note the order and the position of the boundary
x_min and x_max is right


compton error appears in the situation that electron and bw_electron both exist in simulation, leading to a mistake of electron_species

add write_injector_depths, but return if the direction is equal to c_dir_z

some error appears in FRB_1keV: close write_epoch_source_info

add v_inject_z as time and position function

add resonance_scatter but does not benchmark yet
20201013
add resonance_scatter and benchmark and revised a problem in v_inject_z

20201014
some error appears in pair production with photon inject_z
it seems to be a error between electron and bw_electron
20201015
revised!
Epoch原本需要开启produce_photon才能算正负电子对，给breit_wheeler_electron赋值[src/physic.../photons.F90],所以导致bw_electron_id = -1
读取时发生错误，只需要取消掉produce_photon的限制就行
光子的calc_ekbar有问题,问题在于src/housekeeping/partlist/nvar 少了cs_times的项，通讯后发生了错误。

20201017 add positron
20201018 revise relativistic current output problem 
add merge photon of Xu Zheng into 4.17.9
bugs:
    1.对于beam case的处理存在问题：
        1.min_px = max_px 初始设置为是不对的 0.0_num
        2.对于平行的判定利用 == 0.0_num是不行的
    2.计算利用二次方程判断式计算对应的两个粒子的速度解可能会出数值bug，改成矢量的计算。
    3.设置了能量上限，保护能量大于mev的粒子不会被一大块低能光子吸收掉。
    4.准备改成log区间划分

把merge 算法拓展到任意粒子的时候出错
原因在于deallocate报错，可能是由于数组越界
应该是正电子注入时没有给粒子能量给初值导致的
后面发现在正负电子中似乎beam case的合并会带来很多问题，尤其是多删除几个粒子导致能动量出现问题。因此提高了beam判断的要求0.001
还针对极端情况，单能电子入射的情况，判断不合并
merge 算法对正负电子对使用，会导致电磁场不守恒。暂时舍弃。

碰撞模块np在计算e1,e2的时候用的是除了c的p因此能量有问题。需修改。

**********************Polarization****************************
1.一定要开QED模块

**********************Temperature*****************************
1.温度计算里面没有乘weight
