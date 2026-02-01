# SWAT+ vs PWC 农药模型对比分析报告

**项目**: SWAT+与PWC农药计算模型对比研究
**文档版本**: 1.0
**生成日期**: 2026-02-01
**作者**: 基于源代码分析的自动生成报告

---

## 目录

1. [执行摘要](#执行摘要)
2. [模型概述](#模型概述)
3. [核心算法对比](#核心算法对比)
4. [农药施用处理](#农药施用处理)
5. [降解过程模拟](#降解过程模拟)
6. [挥发过程模拟](#挥发过程模拟)
7. [吸附/解吸过程](#吸附解吸过程)
8. [水文过程耦合](#水文过程耦合)
9. [输出结果对比](#输出结果对比)
10. [模型优缺点分析](#模型优缺点分析)

---

## 执行摘要

本报告通过深入分析SWAT+和PWC两个模型的源代码，系统比较了它们在农药计算方面的差异。两个模型在农药运移模拟上有显著不同的设计理念和方法论。

### 主要发现

| 对比维度 | SWAT+ | PWC |
|---------|-------|-----|
| **设计目标** | 流域尺度的综合水文-水质模拟 | 农药风险评估的专用计算工具 |
| **空间尺度** | 分布式HRU/子流域尺度 | 单点/田块尺度 |
| **农药模块** | 集成在综合模型中 | 专门针对农药运移设计 |
| **降解模拟** | 一级动力学，支持代谢产物链 | 支持最多3种化学物质（母体+两代） |
| **挥发模拟** | 简化的水体挥发 | 复杂的冠层阻力模型 |
| **时间步长** | 日步长 | 支持亚日步长细分 |

---

## 模型概述

### SWAT+ 模型概述

**开发机构**: USDA Agricultural Research Service (USDA-ARS) + Texas A&M AgriLife Research

**开源状态**: 完全开源，GitHub仓库: https://github.com/swat-model/swatplus

**核心特性**:
- 基于物理过程的半分布式模型
- 以HRU（水文响应单元）为基本计算单元
- 日时间步长模拟
- 同时模拟水文、泥沙、营养物、农药等多种污染物

**农药模块架构**:
```
pesticide_data_module.f90       - 农药数据库和参数定义
pesticide_init.f90              - 农药模块初始化
ch_pesticide_module.f90         - 河道农药运移模块
hru_pesticide_output.f90        - HRU农药输出
aqu_pesticide_module.f90        - 含水层农药模块
res_pesticide_module.f90        - 水库农药模块
output_ls_pesticide_module.f90  - 农药平衡输出模块
```

### PWC 模型概述

**开发机构**: US EPA (Environmental Protection Agency)

**核心特性**:
- 确定性场景模型
- 基于PRZM（农药根区模型）的田间运移模拟
- 基于VVWM（可变水体水量模型）的水体模拟
- 支持亚日时间步长以提高数值稳定性
- 专门用于监管风险评估

**农药模块架构**:
```
Pesticide_Application.f90       - 农药施用处理
Chemica_Transport.f90           - 化学物质运移核心
Plant_Pesticide_Processes.f90   - 植物农药过程
volatilization.f90              - 挥发过程模拟
vvwm/                           - VVWM水体模块
  ├── CoreCalculations.F90      - 核心计算
  ├── degradation.F90           - 降解过程
  └── VolumeAndWashout.F90      - 水量与冲洗
```

---

## 核心算法对比

### 水文耦合方式

#### SWAT+: 顺序耦合，水文驱动农药运移

SWAT+采用顺序耦合方式，农药运移计算依赖于预先计算的水文变量：

```fortran
! 农药运移使用水文模块输出的变量
use soil_module
use hydrograph_module

! 农药浓度计算基于水文条件
conc = pesticide_mass / (water_volume + soil_bulk_density * Kd)
```

**特点**:
- 水文和农药计算分离
- 计算效率高
- 缺乏双向反馈

#### PWC: 顺序耦合，水文先行 + 溶质隐式求解

PWC水文先行，溶质运移采用隐式求解：

```fortran
! 三对角矩阵求解同时考虑水和农药运移
call setup_tridiagonal_for_przm(subdelt, K, theta_new, theta_old, &
                                 theta_air_new, theta_air_old, old_conc, new_conc)
```

**特点**:
- 先计算水文，再进行溶质运移
- 溶质方程使用隐式三对角/预测-校正器
- 支持亚日时间步长提高数值稳定性

### 数值方法

| 方面 | SWAT+ | PWC |
|------|-------|-----|
| **空间离散** | HRU/土壤层 | 土壤 compartments |
| **时间积分** | 显式欧拉法 | 预测-校正器方法 |
| **对流-扩散** | 简化的一阶近似 | 三对角矩阵求解 |
| **非线性处理** | 线性化处理 | 迭代预测-校正 |

---

## 农药施用处理

### SWAT+ 农药施用

**支持施用方式**（基于源码分析）:
1. 土壤表面施用
2. 叶面施用
3. 土壤混施

**代码示例** ([pesticide_init.f90:67-69](../swatplus/src/pesticide_init.f90#L67-L69)):
```fortran
! 按LAI比例分配叶面农药
pl_frac = pcom(ihru)%plg(ipl)%lai / pcom(ihru)%lai_sum
cs_pl(ihru)%pl_on(ipl)%pest(ipest) = cs_pl(ihru)%pl_on(ipl)%pest(ipest) + &
                                     pl_frac * pest_soil_ini(ipest_db)%plt(ipest)
```

**特点**:
- 简单的施用方式
- 基于LAI分配叶面农药
- 初始浓度直接设定

### PWC 农药施用

**支持施用方式** ([Pesticide_Application.f90:104-216](../PWC_src/Pesticide_Application.f90#L104-L216)):

| CAM代码 | 施用方式 | 描述 |
|---------|----------|------|
| 1 | 土壤表面施用 | 线性递减分布至4cm深度 |
| 2 | 叶面施用 | 部分叶面，部分土壤表面递减分布 |
| 3 | 均匀混施 | 均匀分布至指定深度 |
| 4 | 特定深度施用 | 全部农药置于指定深度 |
| 5 | T-Band施用 | 顶部2cm + 下层均匀分布 |
| 6 | 线性递减施用 | 按深度线性递减 |
| 7 | 线性递增施用 | 按深度线性递增 |
| 8 | 自定义叶面施用 | 用户定义混施深度 |

**关键代码片段** ([Pesticide_Application.f90:326-363](../PWC_src/Pesticide_Application.f90#L326-L363)):
```fortran
SUBROUTINE pesticide_decreasing_distribution(APPAMT,DMAX,BASE,SLOPE,applied_to_soil)
    ! 精确的递减分布算法
    do i=1, ncmpwithchem
        bottom = min(bottom + delx(i), DMAX)
        curr = bottom*bottom*SLOPE/2 + BASE*bottom
        applied_to_soil(i) = (curr - prev) * APPAMT
        prev = curr
    END DO
END SUBROUTINE pesticide_decreasing_distribution
```

**特点**:
- 8种精细的施用方式
- 精确的垂直分布算法
- 支持天气调整施用日期

### 对比总结

| 特性 | SWAT+ | PWC |
|------|-------|-----|
| 施用方式数量 | 3种 | 8种 |
| 垂直分布精度 | 简化处理 | 精确算法 |
| 天气调整 | 无 | 有 (基于降雨窗口) |
| 冠层截留 | 基于LAI | 基于覆盖率 |

---

## 降解过程模拟

### SWAT+ 降解模拟

**降解类型** ([pesticide_data_module.f90:34-40](../swatplus/src/pesticide_data_module.f90#L34-L40)):
```fortran
type pesticide_cp
    real :: decay_f = 0.    ! 叶面降解指数
    real :: decay_s = 0.    ! 土壤降解指数
    real :: decay_a = 0.    ! 水体降解指数
    real :: decay_b = 0.    ! 底泥降解指数
end type pesticide_cp
```

**算法**: 一级动力学
```
C(t) = C0 * exp(-k * t)
```

**输入参数** ([pesticide_data_module.f90:5-22](../swatplus/src/pesticide_data_module.f90#L5-L22)):
- `foliar_hlife`: 叶面半衰期 (天)
- `soil_hlife`: 土壤半衰期 (天)
- `aq_hlife`: 水体半衰期 (天)
- `ben_hlife`: 底泥半衰期 (天)

**代谢产物处理** ([pesticide_data_module.f90:25-32](../swatplus/src/pesticide_data_module.f90#L25-L32)):
```fortran
type daughter_decay_fractions
    real :: foliar_fr = 0.    ! 叶面降解转子产物比例
    real :: soil_fr = 0.      ! 土壤降解转子产物比例
    real :: aq_fr = 0.        ! 水体降解转子产物比例
    real :: ben_fr = 0.       ! 底泥降解转子产物比例
end type daughter_decay_fractions
```

### PWC 降解模拟

**降解类型** ([Chemica_Transport.f90:431-434](../PWC_src/Chemica_Transport.f90#L431-L434)):
```fortran
! 三相降解
DKFLUX(K,I) = concentration(i) * (DWRATE(K,i)*theta_end(i) +       &
                                 DSRATE(K,i)*bulkdensity(i)*kd_new(K,i) + &
                                 DGRATE(K,i)*THAIR_new(i)*new_henry(K,i))
```

**三相降解模型**:
- `DWRATE`: 水相降解速率
- `DSRATE`: 吸附相降解速率
- `DGRATE`: 气相降解速率

**温度校正** ([Chemica_Transport.f90:164-172](../PWC_src/Chemica_Transport.f90#L164-L172)):
```fortran
IF (is_temperature_simulated) THEN
    ! 调整Henry常数
    CALL Henry_Temp_Correction (soil_temp, Henry_unitless(K), ENPY(K), &
                                henry_ref_temp, NCOM2, OLDKH)
    ! Q10温度校正降解速率
    call Q10DK
ENDIF
```

**多级降解产物** ([Plant_Pesticide_Processes.f90:54-103](../PWC_src/Plant_Pesticide_Processes.f90#L54-L103)):
- 支持最多3种化学物质（母体+两代）
- 精确的解析解处理等降解速率情况
- 摩尔转换系数计算

**关键代码** ([Plant_Pesticide_Processes.f90:64-101](../PWC_src/Plant_Pesticide_Processes.f90#L64-L101)):
```fortran
! 二级降解产物
If (nchem == 2 .or. nchem == 3) THEN
    ex2 = EXP((-fol_deg(2))*DELT)
    If (fol_deg(2)==fol_deg(1)) THEN
        r = delt*ex2
    else
        r = (ex1-ex2)/(fol_deg(2)-fol_deg(1))
    END IF
    FOLPST(2) = foliar_pest_initial(1)*foliar_formation_ratio_12*fol_deg(1)*r + &
                foliar_pest_initial(2)*ex2
END IF

! 三级降解产物
If (nchem == 3) THEN
    ! 复杂的解析解计算
    term1 = fol_deg(1)*fol_deg(2)/(fol_deg(2)-fol_deg(1))
    term2 = foliar_formation_ratio_12*foliar_formation_ratio_23*Foliar_Pest_initial(1)
    ...
    FOLPST(3) = term50 + term90 + term100
END IF
```

### 对比总结

| 特性 | SWAT+ | PWC |
|------|-------|-----|
| 降解相数 | 单相（简化） | 三相（水/吸附/气） |
| 温度校正 | 有限 | 完整的Q10校正 |
| 代谢产物 | 任意数量 | 最多3种化学物质（母体+两代） |
| 等速率处理 | 未特别处理 | 精确解析解 |
| Henry常数温度校正 | 无 | 有 |

---

## 挥发过程模拟

### SWAT+ 挥发模拟

**实现方式**: 简化的水体挥发模型

**代码位置**: `aqu_pesticide_module.f90`, `res_pesticide_module.f90`

**参数** ([pesticide_data_module.f90:13](../swatplus/src/pesticide_data_module.f90#L13)):
```fortran
real :: aq_volat = 0.    ! m/day, 水体挥发系数
```

**特点**:
- 仅在水体中考虑挥发
- 使用简化的挥发系数
- 不考虑土壤挥发

### PWC 挥发模拟

**实现方式**: 复杂的冠层阻力模型

**核心算法** ([volatilization.f90:20-131](../PWC_src/volatilization.f90#L20-L131)):

```fortran
subroutine volatilization_setup(k)
    ! 1. 计算空气相扩散系数
    DGAIR = (THAIR_new**(10./3.)/theta_sat**2)*DAIR(K) * THAIR_new

    ! 2. 边界层传导率
    CNDBDY(K) = DAIR(K) / Height_stagnant_air_layer_cm

    ! 3. 如果有冠层，计算冠层阻力
    IF (HEIGHT > Minimum_Canopy_Height_cm) THEN
        Call Get_Crop_Params (zch, z0, d)
        urh = WIND * 864.0
        uch = urh * Log((zch-d)/z0) / Log((uWind_Reference_Height-uWind_D)/uWind_z0)

        IF(Henry_unitless(K).GT.0.0.AND.URH.GT.0.0) THEN
            CALL Canopy(ATEMP, PWIND, ZCH, TOTCR)
            CONDUC(K) = 1.0 / (1.0/CNDBDY(K) + TOTCR)
        END IF
    END IF
end subroutine
```

**冠层阻力模型** ([volatilization.f90:182-350](../PWC_src/volatilization.f90#L182-L350)):
- Richardson数计算大气稳定度
- 稳定度函数修正
- 摩擦速度计算
- 热涡动扩散系数计算
- 冠层分层积分阻力

**挥发通量计算** ([Chemica_Transport.f90:421-423](../PWC_src/Chemica_Transport.f90#L421-L423)):
```fortran
! 土壤表面挥发
PVFLUX(K,1) = -CONDUC(K)*concentration(1)*new_henry(K,1)

! 土壤层间挥发
do i=2,NCOM2-1
    PVFLUX(K,I) = DGAIR(I)*new_henry(K,I)/(0.5*(DELX(I)+DELX(I+1)))*concentration(I) - &
                  DGAIR(I+1)*new_henry(K,I+1)/(0.5*(DELX(I)+DELX(I+1)))*concentration(I+1)
END DO
```

### 对比总结

| 特性 | SWAT+ | PWC |
|------|-------|-----|
| 挥发位置 | 仅水体 | 土壤+水体 |
| 阻力模型 | 简化系数 | 完整冠层阻力模型 |
| 大气稳定度 | 不考虑 | Richardson数计算 |
| 风速剖面 | 不考虑 | 对数风剖面 |
| Henry常数温度校正 | 无 | 有 |
| 空气相扩散 | 不考虑 | 考虑孔隙度和弯曲度 |

---

## 吸附/解吸过程

### SWAT+ 吸附模型

**吸附模型**: 线性吸附（Kd = Koc × OC/100）

**参数** ([pesticide_data_module.f90:7](../swatplus/src/pesticide_data_module.f90#L7)):
```fortran
real :: koc = 0.    ! (mL/g) 土壤有机碳标准化吸附系数
```

**特点**:
- 基于Koc和土壤有机碳含量计算Kd
- 简化的吸附处理
- 不考虑非线性等温吸附

### PWC 吸附模型

**吸附模型**: 支持线性和非线性Freundlich吸附

**线性吸附**:
```
S = Kd * C
```

**非线性Freundlich吸附** ([Chemica_Transport.f90:211-222](../PWC_src/Chemica_Transport.f90#L211-L222)):
```fortran
if (is_Freundlich) THEN
    ! 预测器步骤
    call setup_tridiagonal_for_przm(subdelt, K, theta_new, theta_old, &
                                  theta_air_new, theta_air_old, old_conc, predicted_conc)
    ! 用预测浓度计算新的Kd
    Call Freundlich(k, predicted_conc)
END IF

! 校正器步骤
call setup_tridiagonal_for_przm(subdelt, K, theta_new, theta_old, &
                                 theta_air_new, theta_air_old, old_conc, new_conc)
```

**非平衡吸附** ([NonidealSorption.f90`](../PWC_src/NonidealSorption.f90)):
```fortran
! 双域模型: 移动域和不可移动域
call nonequilibrium(subdelt, k, old_conc, S2_old, theta_old_subday, &
                    theta_air_old_subday, conc1_neq, S2_neq)
```

**特点**:
- 支持预测-校正器方法处理非线性
- 支持非平衡吸附（双域模型）
- 迭代求解保证质量守恒

### 对比总结

| 特性 | SWAT+ | PWC |
|------|-------|-----|
| 吸附模型 | 线性 | 线性 + Freundlich非线性 |
| 非线性处理 | 无 | 预测-校正器迭代 |
| 非平衡吸附 | 无 | 双域模型 |
| 数值稳定性 | 简单处理 | 亚日步长 + 迭代 |

---

## 水文过程耦合

### SWAT+ 水文耦合

**耦合方式**: 顺序耦合

**流程**:
```
水文模块 → 径流/渗流/蒸散发计算 → 农药运移使用水文结果
```

**特点**:
- 计算效率高
- 模块化设计
- 缺乏双向反馈

### PWC 水文耦合

**耦合方式**: 顺序耦合（先水文后溶质）

**三对角矩阵系统** ([Chemica_Transport.f90:266-382](../PWC_src/Chemica_Transport.f90#L266-L382)):
```fortran
SUBROUTINE setup_tridiagonal_for_przm(delt, K, theta_new, theta_old, ...)
    ! a(i)y(i-1) + b(i)y(i) + c(i)y(i+1) = f(i)

    ! 对流项
    B(I) = ... + VEL(I)*theta_new(I)/DELX(I) + ...

    ! 扩散项
    B(I) = ... + (dispersion(I)*theta_new(I) + new_henry(K,I)*DGAIR(I))/DELX(I)**2 + ...

    ! 降解项
    B(I) = ... + (DWRATE(K,I)*theta_new(I)) + ...
                (DSRATE(K,I)*kd_new(K,I)*bulkdensity(I)) + ...
                (DGRATE(K,I)*theta_air_new(I)*new_henry(K,I)) + ...

    ! 径流项
    B(I) = ... + runoff_on_day*runoff_intensity(i) + ...

    ! 侵蚀项
    B(I) = ... + enriched_eroded_solids*kd_new(K,i)*erosion_intensity(i) + ...

    ! 求解三对角系统
    CALL TRIDIAGONAL_Solution (A,B,C, new_conc,F,NCOM2)
END SUBROUTINE
```

**亚日时间步长** ([Chemica_Transport.f90:184-236](../PWC_src/Chemica_Transport.f90#L184-L236)):
```fortran
do j = 1, number_subdelt
    ! 水分和空气含量在亚日步长内变化
    theta_new_subday = theta_old_subday + delta_watercontent
    theta_air_new_subday = theta_air_old_subday + delta_aircontent

    ! 预测-校正器求解
    call setup_tridiagonal_for_przm(subdelt, K, ...)
END DO
```

### 对比总结

| 特性 | SWAT+ | PWC |
|------|-------|-----|
| 耦合方式 | 顺序 | 顺序（溶质隐式求解） |
| 时间步长 | 日步长 | 支持亚日细分 |
| 数值方法 | 简化近似 | 三对角矩阵求解 |
| 双向反馈 | 无 | 无（溶质不反作用于水文） |
| 计算效率 | 高 | 较低 |
| 数值精度 | 较低 | 较高 |

---

## 输出结果对比

### SWAT+ 输出变量

**HRU农药平衡** ([output_ls_pesticide_module.f90:5-24](../swatplus/src/output_ls_pesticide_module.f90#L5-L24)):
```fortran
type pesticide_balance
    real :: plant = 0.      ! 植物上农药 (kg/ha)
    real :: soil = 0.       ! 土壤中农药 (kg/ha)
    real :: sed = 0.        ! 随泥沙流失 (kg/ha)
    real :: surq = 0.       ! 地表径流流失 (kg/ha)
    real :: latq = 0.       ! 侧向流流失 (kg/ha)
    real :: tileq = 0.      ! 排水流失 (kg/ha)
    real :: perc = 0.       ! 淋溶流失 (kg/ha)
    real :: apply_s = 0.    ! 土壤施用 (kg/ha)
    real :: apply_f = 0.    ! 叶面施用 (kg/ha)
    real :: decay_s = 0.    ! 土壤降解 (kg/ha)
    real :: decay_f = 0.    ! 叶面降解 (kg/ha)
    real :: wash = 0.       ! 冲刷 (kg/ha)
    real :: metab_s = 0.    ! 土壤代谢 (kg/ha)
    real :: metab_f = 0.    ! 叶面代谢 (kg/ha)
    real :: pl_uptake = 0.  ! 植物吸收 (kg/ha)
    real :: in_plant = 0.   ! 植物内农药 (kg/ha)
end type pesticide_balance
```

**河道农药过程** ([ch_pesticide_module.f90:10-25](../swatplus/src/ch_pesticide_module.f90#L10-L25)):
```fortran
type ch_pesticide_processes
    real :: tot_in = 0.       ! 总入流 (kg)
    real :: sol_out = 0.      ! 溶解态出流 (kg)
    real :: sor_out = 0.      ! 吸附态出流 (kg)
    real :: react = 0.        ! 反应损失 (kg)
    real :: metab = 0.        ! 代谢 (kg)
    real :: volat = 0.        ! 挥发 (kg)
    real :: settle = 0.       ! 沉降 (kg)
    real :: resus = 0.        ! 再悬浮 (kg)
    real :: difus = 0.        ! 扩散 (kg)
    real :: react_bot = 0.    ! 底泥反应 (kg)
    real :: metab_bot = 0.    ! 底泥代谢 (kg)
    real :: bury = 0.         ! 底泥掩埋 (kg)
    real :: water = 0.        ! 水体农药量 (kg)
    real :: benthic = 0.      ! 底泥农药量 (kg)
end type ch_pesticide_processes
```

### PWC 输出变量

**田间通量** ([Chemica_Transport.f90:386-488](../PWC_src/Chemica_Transport.f90#L386-L488)):
```fortran
! 主要通量
ROFLUX(K)     ! 径流通量 (g/cm²/day)
ERFLUX(K)     ! 侵蚀通量 (g/cm²/day)
PVFLUX(K,:)   ! 挥发通量 (g/cm²/day)
DKFLUX(K,:)   ! 降解通量 (g/cm²/day)
WOFLUX(K)     ! 植物冲刷通量 (g/cm²/day)
UPFLUX(K,:)   ! 植物吸收通量 (g/cm²/day)
DCOFLUX(K)    ! 底部出流通量 (g/cm²/day)
SRCFLX(K,:)   ! 降解产物生成通量 (g/cm²/day)

! 详细分解
SDKFLX(K)     ! 总降解通量
SUPFLX(K)     ! 总植物吸收通量
```

### 对比总结

| 特性 | SWAT+ | PWC |
|------|-------|-----|
| 输出粒度 | HRU/河道/子流域 | 田间单点 |
| 空间输出 | 分布式 | 单点 |
| 时间分辨率 | 日/月/年 | 日 |
| 质量平衡 | 完整 | 完整 |
| 通量分解 | 中等 | 详细 |
| 中间变量 | 有限 | 丰富 |

---

## 模型优缺点分析

### SWAT+ 优点

1. **综合性强**: 同时模拟多种污染物和过程
2. **分布式**: 支持流域尺度的空间分布
3. **开源**: 完全开源，可自由修改和扩展
4. **社区支持**: 庞大的用户群和丰富的文献
5. **计算效率**: 日步长和简化算法保证效率
6. **数据需求**: 相对较低，适合数据稀缺地区

### SWAT+ 缺点

1. **农药模块简化**: 挥发、吸附等过程简化处理
2. **时间步长限制**: 固定日步长，无法捕捉短时事件
3. **双向耦合缺失**: 农药不影响水文过程
4. **数值精度**: 简化算法可能损失精度
5. **校准需求**: 参数众多，校准工作量大

### PWC 优点

1. **专业性强**: 专门针对农药风险评估设计
2. **算法精细**: 完整的物理过程模拟
3. **数值稳定**: 预测-校正器和亚日步长保证稳定性
4. **监管认可**: EPA官方模型，用于监管决策
5. **过程完整**: 考虑三相、非线性、非平衡等复杂过程
6. **文档完善**: 详细的技术文档和用户手册

### PWC 缺点

1. **单点模型**: 不支持空间分布
2. **非完全开源**: 软件需官方获取，源码不完整公开
3. **计算成本**: 精细算法导致计算时间较长
4. **学习曲线**: 复杂的理论基础，学习难度大
5. **数据需求**: 需要详细的输入参数
6. **扩展性差**: 难以与其他模型集成

---

## 应用建议

### 选择SWAT+的场景

1. **流域尺度研究**: 需要空间分布结果
2. **综合评估**: 同时考虑多种污染物
3. **数据有限**: 参数数据不完整
4. **长期模拟**: 多年连续模拟
5. **情景分析**: 多种管理措施比较
6. **开源需求**: 需要修改或扩展模型

### 选择PWC的场景

1. **风险评估**: 农药注册和再评估
2. **监管合规**: 需要EPA认可的方法
3. **精度要求**: 高精度的浓度预测
4. **单点详细**: 田间尺度的精细模拟
5. **短期事件**: 暴雨径流等短时事件
6. **复杂化学**: 需要考虑非线性、非平衡过程

---

## 文献来源

1. [Extending SWAT2012 and SWAT+ models (HESS, 2025)](https://hess.copernicus.org/articles/29/6703/2025/)
2. [PWC User Manual (EPA, 2015)](https://19january2021snapshot.epa.gov/sites/static/files/2015-12/documents/pwc_user_manual_12-8-15.pdf)
3. [Advancing Pesticide Exposure Assessments with SWAT+](https://swat.tamu.edu/media/t4ijyo52/a3-1-swatplus_strassbourg_03jul2024.pdf)
4. [PWC tool for soil concentration estimates](https://academic.oup.com/ieam/article/21/4/943/8043221)
5. [SWAT+ Official Documentation](https://swatplus.gitbook.io/docs)

---

*报告结束*
