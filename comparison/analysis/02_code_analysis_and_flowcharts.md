# SWAT+ vs PWC 农药模块代码分析与流程图

**文档版本**: 1.0
**生成日期**: 2026-02-01

---

## 目录

1. [SWAT+ 农药模块代码分析](#swat-农药模块代码分析)
2. [PWC 农药模块代码分析](#pwc-农药模块代码分析)
3. [算法流程图对比](#算法流程图对比)
4. [关键函数对比](#关键函数对比)

---

## SWAT+ 农药模块代码分析

### 模块架构图

```
┌─────────────────────────────────────────────────────────────────┐
│                        SWAT+ 主程序                              │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                      pesticide_init                             │
│  - 分配HRU农药数组                                               │
│  - 初始化土壤和植物农药浓度                                       │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                   每日时间步长循环                                │
└─────────────────────────────────────────────────────────────────┘
                              │
              ┌───────────────┼───────────────┐
              ▼               ▼               ▼
     ┌─────────────┐  ┌─────────────┐  ┌─────────────┐
     │ HRU 农药计算 │  │ 河道农药运移│  │ 水库农药模拟│
     └─────────────┘  └─────────────┘  └─────────────┘
              │               │               │
              ▼               ▼               ▼
     ┌─────────────────────────────────────────────────┐
     │            pesticide_balance 输出                │
     │  - plant(植物)                                   │
     │  - soil(土壤)                                    │
     │  - surq(地表径流)                                │
     │  - latq(侧向流)                                  │
     │  - perc(淋溶)                                    │
     │  - decay_s/decay_f(降解)                         │
     └─────────────────────────────────────────────────┘
```

### 核心数据结构分析

#### 1. 农药数据库类型 ([pesticide_data_module.f90:5-22](../swatplus/src/pesticide_data_module.f90#L5-L22))

```fortran
type pesticide_db
    character(len=16) :: name = ""      ! 农药名称
    real :: koc = 0.                    ! 土壤有机碳吸附系数 (mL/g)
    real :: washoff = 0.                ! 叶面冲刷系数
    real :: foliar_hlife = 0.           ! 叶面半衰期 (天)
    real :: soil_hlife = 0.             ! 土壤半衰期 (天)
    real :: solub = 0.                  ! 水溶性 (mg/L)
    real :: aq_hlife = 0.               ! 水体半衰期 (天)
    real :: aq_volat = 0.               ! 水体挥发系数 (m/day)
    real :: mol_wt = 0.                 ! 分子量 (g/mol)
    real :: aq_resus = 0.               ! 水体再悬浮速度 (m/day)
    real :: aq_settle = 0.              ! 水体沉降速度 (m/day)
    real :: ben_act_dep = 0.            ! 底泥活跃层深度 (m)
    real :: ben_bury = 0.               ! 底泥掩埋速度 (m/day)
    real :: ben_hlife = 0.              ! 底泥半衰期 (天)
    real :: pl_uptake = 0.              ! 植物吸收比例
end type pesticide_db
```

**关键参数说明**:
- `koc`: 土壤吸附的关键参数，通过 `Kd = Koc * organic_carbon` 计算分配系数
- `washoff`: 降雨冲刷系数，0-1之间，决定叶面农药被冲刷的比例
- `foliar_hlife/soil_hlife`: 一级降解动力学参数，降解常数 `k = ln(2)/half_life`
- `aq_volat`: 简化的水体挥发系数

#### 2. 农药平衡类型 ([output_ls_pesticide_module.f90:5-24](../swatplus/src/output_ls_pesticide_module.f90#L5-L24))

```fortran
type pesticide_balance
    real :: plant = 0.      ! 植物表面农药量 (kg/ha)
    real :: soil = 0.       ! 土壤中农药总量 (kg/ha)
    real :: sed = 0.        ! 随泥沙流失的农药 (kg/ha)
    real :: surq = 0.       ! 地表径流流失 (kg/ha)
    real :: latq = 0.       ! 侧向流流失 (kg/ha)
    real :: tileq = 0.      ! 排水流失 (kg/ha)
    real :: perc = 0.       ! 淋溶流失 (kg/ha)
    real :: apply_s = 0.    ! 土壤施用量 (kg/ha)
    real :: apply_f = 0.    ! 叶面施用量 (kg/ha)
    real :: decay_s = 0.    ! 土壤降解量 (kg/ha)
    real :: decay_f = 0.    ! 叶面降解量 (kg/ha)
    real :: wash = 0.       ! 叶面冲刷量 (kg/ha)
    real :: metab_s = 0.    ! 土壤代谢产物量 (kg/ha)
    real :: metab_f = 0.    ! 叶面代谢产物量 (kg/ha)
    real :: pl_uptake = 0.  ! 植物吸收量 (kg/ha)
    real :: in_plant = 0.   ! 植物体内农药量 (kg/ha)
end type pesticide_balance
```

### 农药初始化流程

```fortran
! 文件: pesticide_init.f90
subroutine pesticide_init
    ! 1. 分配HRU农药数组
    do ihru = 1, sp_ob%hru
        nly = soil(ihru)%nly              ! 土壤层数
        npl = pcom(ihru)%npl              ! 植物数量

        ! 分配土壤层农药
        allocate (cs_soil(ihru)%ly(nly))
        do ly = 1, nly
            allocate (cs_soil(ihru)%ly(ly)%pest(npmx))
        end do

        ! 分配植物农药
        allocate (cs_pl(ihru)%pl_in(npl))   ! 植物内农药
        allocate (cs_pl(ihru)%pl_on(npl))   ! 植物表面农药
        allocate (cs_pl(ihru)%pl_up(npl))   ! 植物吸收
    end do

    ! 2. 设置初始农药浓度
    do ipest = 1, npmx
        ! 植物表面初始农药 (按LAI分配)
        do ipl = 1, pcom(ihru)%npl
            pl_frac = pcom(ihru)%plg(ipl)%lai / pcom(ihru)%lai_sum
            cs_pl(ihru)%pl_on(ipl)%pest(ipest) = pl_frac * initial_amount
        end do

        ! 土壤初始农药 (均匀分布各层)
        do ly = 1, soil(ihru)%nly
            wt1 = soil(ihru)%phys(ly)%bd * soil(ihru)%phys(ly)%thick / 100.
            cs_soil(ihru)%ly(ly)%pest(ipest) = solpst * wt1
        end do
    end do
end subroutine
```

### HRU农药输出流程

```fortran
! 文件: hru_pesticide_output.f90
subroutine hru_pesticide_output(ihru)
    do ipest = 1, cs_db%num_pests
        ! 累积月度和年度数据
        hpestb_m(j)%pest(ipest) = hpestb_m(j)%pest(ipest) + hpestb_d(j)%pest(ipest)

        ! 日输出
        if (pco%day_print == "y") then
            write (2800,100) time%day, time%mo, time%day_mo, time%yrc, &
                           j, ob(iob)%gis_id, ob(iob)%name, &
                           cs_db%pests(ipest), hpestb_d(j)%pest(ipest)
        end if

        ! 月输出
        if (time%end_mo == 1) then
            hpestb_y(j)%pest(ipest) = hpestb_y(j)%pest(ipest) + hpestb_m(j)%pest(ipest)
            const = float (ndays(time%mo + 1) - ndays(time%mo))
            hpestb_m(j)%pest(ipest) = hpestb_m(j)%pest(ipest) // const
        end if

        ! 年输出
        if (time%end_yr == 1) then
            hpestb_a(j)%pest(ipest) = hpestb_a(j)%pest(ipest) + hpestb_y(j)%pest(ipest)
        end if
    end do
end subroutine
```

### 河道农药过程

```fortran
! 文件: ch_pesticide_module.f90
type ch_pesticide_processes
    real :: tot_in = 0.       ! 总入流量 (kg)
    real :: sol_out = 0.      ! 溶解态出流 (kg)
    real :: sor_out = 0.      ! 吸附态出流 (kg)
    real :: react = 0.        ! 反应损失 (kg)
    real :: metab = 0.        ! 代谢量 (kg)
    real :: volat = 0.        ! 挥发量 (kg)
    real :: settle = 0.       ! 沉降量 (kg)
    real :: resus = 0.        ! 再悬浮量 (kg)
    real :: difus = 0.        ! 扩散量 (kg)
    real :: react_bot = 0.    ! 底泥反应 (kg)
    real :: metab_bot = 0.    ! 底泥代谢 (kg)
    real :: bury = 0.         ! 底泥掩埋 (kg)
    real :: water = 0.        ! 水体农药量 (kg)
    real :: benthic = 0.      ! 底泥农药量 (kg)
end type ch_pesticide_processes
```

---

## PWC 农药模块代码分析

### 模块架构图

```
┌─────────────────────────────────────────────────────────────────┐
│                          PWC 主程序                              │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                       Initialization                             │
│  - ReadInputs.f90: 读取输入文件                                  │
│  - 调整施用日期 (天气窗口)                                        │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                    日循环 (chem_transport_onfield)               │
└─────────────────────────────────────────────────────────────────┘
                              │
              ┌───────────────┼───────────────┐
              ▼               ▼               ▼
     ┌─────────────┐  ┌─────────────┐  ┌─────────────┐
     │ 农药施用     │  │ 植物农药处理 │  │ 田间运移    │
     │ PESTAP      │  │ Plant_      │  │ Chemical_   │
     │             │  │ Pesticide_  │  │ Transport   │
     │             │  │ Processes   │  │             │
     └─────────────┘  └─────────────┘  └─────────────┘
                               │
                               ▼
     ┌─────────────────────────────────────────────────┐
     │          Chemical_Transport_Oneday              │
     │  1. 检查施用日期 → PESTAP                       │
     │  2. 植物农药冲刷                                │
     │  3. 植物农药降解                                │
     │  4. 田间土壤运移 (亚日步长循环)                  │
     │     a. 非平衡吸附                                │
     │     b. 预测器 (Freundlich)                      │
     │     c. 校正器                                    │
     │  5. 通量计算                                    │
     │  6. 输出到VVWM                                  │
     └─────────────────────────────────────────────────┘
                               │
                               ▼
     ┌─────────────────────────────────────────────────┐
     │                   VVWM水体模块                   │
     │  - 水体农药浓度计算                              │
     │  - 降解、挥发、沉降、再悬浮                      │
     │  - 底泥-水交换                                  │
     └─────────────────────────────────────────────────┘
```

### 核心数据结构分析

#### 1. 农药施用方式 ([Pesticide_Application.f90:104-216](../PWC_src/Pesticide_Application.f90#L104-L216))

PWC支持8种精细的农药施用方式:

| CAM | 方式 | 描述 | 代码位置 |
|-----|------|------|----------|
| 1 | 土壤表面 | 线性递减至4cm | case 1 |
| 2 | 叶面施用 | 部分叶面，部分递减 | case 2 |
| 3 | 均匀混施 | 均匀分布至DEPI深度 | case 3 |
| 4 | 特定深度 | 全部置于DEPI深度 | case 4 |
| 5 | T-Band | 顶部2cm + 下层均匀 | case 5 |
| 6 | 递减施用 | 深度线性递减 | case 6 |
| 7 | 递增施用 | 深度线性递增 | case 7 |
| 8 | 自定义叶面 | 可定义混施深度 | case 8 |

#### 2. 垂直分布算法

**线性递减分布** ([Pesticide_Application.f90:326-363](../PWC_src/Pesticide_Application.f90#L326-L363)):

```fortran
SUBROUTINE pesticide_decreasing_distribution(APPAMT, DMAX, BASE, SLOPE, applied_to_soil)
    ! 输入:
    !   APPAMT  - 总施用量 (g/cm²)
    !   DMAX    - 分布深度 (cm)
    !   BASE    - 初始分布系数
    !   SLOPE   - 递减斜率

    ! 算法: 三角形分布，面积归一化
    !       f(z) = SLOPE*z + BASE
    !       ∫f(z)dz = 1 (归一化)

    DEP = 0.0
    applied_to_soil = 0.0
    prev = 0.0
    bottom = 0.0

    ! 1. 计算受影响的层数
    do ncmpwithchem = 1, ncom2
        DEP = DEP + delx(ncmpwithchem)
        if (DEP >= (DMAX-1E-7)) exit
    END DO

    ! 2. 精确积分分配
    do i=1, ncmpwithchem
        bottom = min(bottom + delx(i), DMAX)
        curr = bottom*bottom*SLOPE/2 + BASE*bottom  ! 积分
        applied_to_soil(i) = (curr - prev) * APPAMT   ! 差值分配
        prev = curr
    END DO
END SUBROUTINE
```

**T-Band分布** ([Pesticide_Application.f90:543-584](../PWC_src/Pesticide_Application.f90#L543-L584)):

```fortran
SUBROUTINE pesticide_Tband_distribution(APPAMT, DMAX, top2cm, applied_to_soil)
    ! T-Band: 顶部2cm集中 + 下层均匀分布
    ! top2cm: 顶部2cm的比例 (0-1)

    do i=1, ncmpwithchem
        bottom = min(bottom + delx(i), DMAX)
        if (bottom <= 2.0) THEN
            ! 顶部2cm: 线性分布
            curr = (bottom/2.0)*APPAMT * top2cm
            applied_to_soil(i) = (curr - prev)
        else
            ! 下层: 均匀分布剩余部分
            curr = APPAMT * top2cm + ((bottom-2.0)/(DMAX-2.0)) * APPAMT * (1-top2cm)
            applied_to_soil(i) = (curr - prev)
        END IF
        prev = curr
    END DO
END SUBROUTINE
```

### 田间化学运移核心算法

**主循环结构** ([Chemica_Transport.f90:84-263](../PWC_src/Chemica_Transport.f90#L84-L263)):

```fortran
SUBROUTINE chemical_transport_oneday
    ! 1. 农药施用检查
    do i = 1, total_applications
        if (application_date(i)==julday1900) THEN
            CALL PESTAP(i)  ! 应用农药
        END IF
    END DO

    ! 2. 植物农药过程
    if (some_applications_were_foliar) THEN
        CALL plant_pesticide_washoff      ! 冲刷
        if (harvest_day) CALL plant_pesticide_harvest_application
        CALL plant_pesticide_degradation  ! 降解
    END IF

    ! 3. 每个化学物质的运移
    do K=1, NCHEM
        ! 3.1 计算初始孔隙水浓度
        do concurrent (I=1:NCOM2)
            old_conc(i) = conc_total_per_water(k,i) * theta_zero(i) / &
                         (theta_zero(i) + kd_new(k,i)*bulkdensity(i) + &
                          thair_old(i)*old_henry(k,i))
        END DO

        ! 3.2 温度校正
        IF (is_temperature_simulated) THEN
            CALL Henry_Temp_Correction(...)
            call Q10DK  ! 降解速率温度校正
        ENDIF

        ! 3.3 挥发设置
        call volatilization_setup(k)

        ! 3.4 亚日时间步长循环
        do j = 1, number_subdelt
            ! 非平衡吸附
            if (is_nonequilibrium) THEN
                call nonequilibrium(subdelt, k, old_conc, S2_old, ...)
            END IF

            ! 预测器 (Freundlich非线性)
            if (is_Freundlich) THEN
                call setup_tridiagonal_for_przm(subdelt, K, ..., predicted_conc)
                Call Freundlich(k, predicted_conc)  ! 更新Kd
            END IF

            ! 校正器
            call setup_tridiagonal_for_przm(subdelt, K, ..., new_conc)

            ! 更新状态
            theta_old_subday = theta_new_subday
            theta_air_old_subday = theta_air_new_subday
        END DO

        ! 3.5 存储结果
        do i=1, NCOM2
            conc_total_per_water(K,I) = new_conc(i) * (theta_end(i) + ...) / theta_end(i)
            mass_in_compartment(K,I) = new_conc(i) * (...) * delx(i)
        END DO

        ! 3.6 通量计算
        call flux_calculations(K, new_conc)
    END DO

    ! 4. 输出到VVWM
    call get_inputs_from_field_for_vvwm
    CALL WRITE_outputfile_2
END SUBROUTINE
```

### 三对角矩阵求解

**矩阵系统** ([Chemica_Transport.f90:266-382](../PWC_src/Chemica_Transport.f90#L266-L382)):

```fortran
SUBROUTINE setup_tridiagonal_for_przm(delt, K, theta_new, theta_old, ...)
    ! 求解系统: a(i)y(i-1) + b(i)y(i) + c(i)y(i+1) = f(i)

    ! 表面层系数
    A(1) = 0.0
    B(1) = ( dispersion*theta/DX²                    ! 扩散
           + VEL*theta/DX                            ! 对流
           + DWRATE*theta                            ! 水相降解
           + DSRATE*Kd*bd                            ! 吸附相降解
           + DGRATE*theta_air*Henry                  ! 气相降解
           + GAMMA1*ET*theta/SoilWater               ! 植物吸收
           + runoff*runoff_intensity                 ! 径流
           + enriched_eroded*Kd*erosion_intensity)   ! 侵蚀
           * DELT
           + theta + Kd*bd + theta_air*Henry         ! 存储项
           + CONDUC*Henry*DELT/DX                    ! 挥发

    C(1) = -(dispersion*theta + Henry*DGAIR) * DELT / (DX*0.5*(DX+DX_next))

    F(1) = (theta_old + Kd_old*bd + theta_air_old*Henry_old) * old_conc
           + washoff*DELT/DX                         ! 冲刷源项

    ! 中间层系数
    do I=2, NCOM2-1
        A(I) = ( -(dispersion_prev*theta_prev + Henry_prev*DGAIR_prev) / ... ) * DELT

        B(I) = ( (dispersion*theta + Henry*DGAIR) / (DX*0.5*(DX_prev+DX))
               + (dispersion*theta + Henry*DGAIR) / (DX*0.5*(DX+DX_next))
               + VEL*theta/DX
               + DWRATE*theta
               + DSRATE*Kd*bd
               + DGRATE*theta_air*Henry
               + GAMMA1*ET*theta/SoilWater
               + runoff*runoff_intensity
               + enriched_eroded*Kd*erosion_intensity) * DELT
               + theta + Kd*bd + theta_air*Henry

        C(I) = -(dispersion_next*theta_next + Henry_next*DGAIR_next) * DELT / ...

        F(I) = (theta_old + Kd_old*bd + theta_air_old*Henry_old) * old_conc
               + washoff*DELT/DX
    END DO

    ! 底层系数 (无径流和侵蚀)
    A(NCOM2) = ...
    B(NCOM2) = ...
    C(NCOM2) = 0.0
    F(NCOM2) = ...

    ! 求解三对角系统
    CALL TRIDIAGONAL_Solution (A, B, C, new_conc, F, NCOM2)
END SUBROUTINE
```

**物理意义**:
- `A(i)`: 上层对当前层的影响
- `B(i)`: 当前层的系数（存储+各种损失项）
- `C(i)`: 下层对当前层的影响
- `F(i)`: 源项（上一步浓度+冲刷输入）

### 通量计算

**各类通量** ([Chemica_Transport.f90:386-488](../PWC_src/Chemica_Transport.f90#L386-L488)):

```fortran
subroutine flux_calculations(k, concentration)
    ! 1. 径流通量
    ROFLUX(K) = sum( runoff_on_day * runoff_intensity(1:RNCMPT) *
                    concentration(1:RNCMPT) * DELX(1:RNCMPT) )

    ! 2. 侵蚀通量
    ERFLUX(K) = sum( enriched_eroded_solids * kd_new(K,1:erosion_compt) *
                    concentration(1:erosion_compt) *
                    erosion_intensity(1:erosion_compt) * DELX(1:erosion_compt) )

    ! 3. 挥发通量
    PVFLUX(K,1) = -CONDUC(K) * concentration(1) * new_henry(K,1)
    do i=2, NCOM2-1
        PVFLUX(K,i) = DGAIR(i)*new_henry(K,i)/(0.5*(DELX(i)+DELX(i+1))) * concentration(i) -
                      DGAIR(i+1)*new_henry(K,i+1)/(0.5*(DELX(i)+DELX(i+1))) * concentration(i+1)
    END DO

    ! 4. 降解通量
    do i=1, NCOM2
        DKFLUX(K,i) = DELX(i) * concentration(i) *
                      ( DWRATE(K,i) * theta_end(i) +
                        DSRATE(K,i) * bulkdensity(i) * kd_new(K,i) +
                        DGRATE(K,i) * THAIR_new(i) * new_henry(K,i) )
    END DO

    ! 5. 植物吸收通量
    UPFLUX(K,1:ncom2) = GAMMA1(K,1:ncom2) * EvapoTran(1:ncom2) * concentration(1:ncom2)

    ! 6. 底部出流通量
    DCOFLUX(K) = VEL(NCOM2) * concentration(NCOM2) * theta_end(NCOM2)

    ! 7. 降解产物生成通量
    if (nchem > 1) then
        if (k == 1) then  ! 母化合物
            srcflx(2,i) = DELX(i) * concentration(i) *
                          ( MolarConvert_aq12(i) * DWRATE(K,i) * theta_end(i) +
                            MolarConvert_s12(i) * DSRATE(K,i) * bulkdensity(i) * kd_new(K,i) )
        end if
    end if
end subroutine
```

---

## 算法流程图对比

### 农药施用流程对比

```
                    SWAT+                                   PWC
                      │                                      │
                      ▼                                      ▼
             ┌───────────────┐                    ┌─────────────────────┐
             │   读取施用    │                    │  读取施用参数       │
             │   操作数据    │                    │  (CAM, DEPI等)     │
             └───────────────┘                    └─────────────────────┘
                      │                                      │
                      ▼                                      ▼
             ┌───────────────┐                    ┌─────────────────────┐
             │ 按LAI分配     │                    │ 天气调整施用日期    │
             │ 叶面农药      │                    │ (雨量窗口检查)      │
             └───────────────┘                    └─────────────────────┘
                      │                                      │
                      ▼                                      ▼
             ┌───────────────┐                    ┌─────────────────────┐
             │ 土壤均匀      │                    │ 根据CAM选择分布     │
             │ 分布          │                    │ 算法               │
             └───────────────┘                    └─────────────────────┘
                      │                                      │
                      ▼                                      ▼
             ┌───────────────┐                    ┌─────────────────────┐
             │ 设置初始浓度  │                    │ 精确积分计算各层    │
             │              │                    │ 农药分布量          │
             └───────────────┘                    └─────────────────────┘
```

### 降解计算流程对比

```
                    SWAT+                                   PWC
                      │                                      │
                      ▼                                      ▼
             ┌───────────────┐                    ┌─────────────────────┐
             │ 读取半衰期    │                    │ 读取三相降解速率    │
             │ (hlife)       │                    │ (水/吸附/气)        │
             └───────────────┘                    └─────────────────────┘
                      │                                      │
                      ▼                                      ▼
             ┌───────────────┐                    │ ┌───────────────────┐ │
             │ k = ln2/hlife │                    │ │ 温度校正?          │ │
             └───────────────┘                    │ └───────────────────┘ │
                      │                            │          │          │
                      ▼                            │    ┌─────┴─────┐    │
             ┌───────────────┐                    │    ▼           ▼    │
             │ C = C0*exp(-k │                    │ │ 是          否   │ │
             │    *dt)       │                    │ ▼             │    │
             └───────────────┘                    │ ┌─────────┐   │    │
                      │                            │ │ Q10     │   │    │
                      ▼                            │ │ 校正    │   │    │
             ┌───────────────┐                    │ └─────────┘   │    │
             │ 简单一级      │                    │      │         │    │
             │ 动力学        │                    │      ▼         ▼    │
             └───────────────┘                    │ ┌───────────────────┐ │
                                                           ▼            │
                                                  ┌─────────────────────┐ │
                                                  │ 解析解计算          │ │
                                                  │ (支持等速率处理)    │ │
                                                  └─────────────────────┘ │
                                                                   │    │
                                                                   ▼    ▼
                                                  ┌──────────────────────┐
                                                  │ 计算母体和降解产物   │
                                                  │ 浓度变化             │
                                                  └──────────────────────┘
```

### 挥发计算流程对比

```
                    SWAT+                                   PWC
                      │                                      │
                      ▼                                      ▼
             ┌───────────────┐                    ┌─────────────────────┐
             │ 读取挥发系数  │                    │ 读取Henry常数       │
             │ (aq_volat)    │                    │ 和空气扩散系数      │
             └───────────────┘                    └─────────────────────┘
                      │                                      │
                      ▼                                      ▼
             ┌───────────────┐                    ┌─────────────────────┐
             │ 仅计算水体    │                    │ 温度校正Henry常数   │
             │ 挥发          │                    │ (参考温度校正)      │
             └───────────────┘                    └─────────────────────┘
                      │                                      │
                      ▼                                      ▼
             ┌───────────────┐                    │ ┌───────────────────┐ │
             │ Flux = -k * C │                    │ │ 有冠层?           │ │
             └───────────────┘                    │ └───────────────────┘ │
                                                      │          │
                                                      ▼          ▼
                                        ┌─────────────────┐ ┌─────────────────┐
                                        │ 是: 冠层阻力   │ │ 否: 边界层阻力  │
                                        │ 模型           │ │ 模型            │
                                        └─────────────────┘ └─────────────────┘
                                                      │          │
                                                      └─────┬─────┘
                                                            ▼
                                        ┌─────────────────────────────┐
                                        │ 计算边界层传导率             │
                                        │ CNDBDY = DAIR / d           │
                                        └─────────────────────────────┘
                                                            │
                                                            ▼
                                        ┌─────────────────────────────┐
                                        │ 计算总传导率                 │
                                        │ CONDUC = 1/(1/CNDBDY + R_can)│
                                        └─────────────────────────────┘
                                                            │
                                                            ▼
                                        ┌─────────────────────────────┐
                                        │ Flux = -CONDUC * C * H      │
                                        └─────────────────────────────┘
```

---

## 关键函数对比

### 农药施用函数

| 特性 | SWAT+ | PWC |
|------|-------|-----|
| **函数名** | `pesticide_init` | `PESTAP` |
| **施用方式** | 3种 | 8种 |
| **垂直分布** | 简化均匀 | 精确积分 |
| **冠层截留** | 基于LAI | 基于覆盖率 |
| **天气调整** | 无 | 有 |

### 降解计算函数

| 特性 | SWAT+ | PWC |
|------|-------|-----|
| **降解模型** | 单相一级动力学 | 三相一级动力学 |
| **温度校正** | 有限 | 完整Q10校正 |
| **等速率处理** | 无 | 精确解析解 |
| **代谢产物** | 简单比例 | 精确摩尔转换 |

### 挥发计算函数

| 特性 | SWAT+ | PWC |
|------|-------|-----|
| **函数名** | 内嵌计算 | `volatilization_setup` |
| **阻力模型** | 简化系数 | 冠层+边界层 |
| **大气稳定度** | 无 | Richardson数 |
| **Henry校正** | 无 | 温度校正 |

---

*文档结束*
