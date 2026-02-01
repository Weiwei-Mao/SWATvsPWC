# SWAT+ vs PWC 农药模型对比项目规划

## 项目概述

**最终目标**: 基于源代码分析，制作一份对比 SWAT+ 和 PWC 农药运移模型的演示文稿（PPT）.PPT使用Warp创建。

本项目旨在对比分析两个农业农药运移计算模型：
- **SWAT+** (Soil and Water Assessment Tool Plus) - USDA-ARS和Texas A&M AgriLife Research联合开发的开源流域尺度模型
- **PWC** (Pesticide in Water Calculator) - EPA开发的农药风险评估模型

## 项目目标

对比两个模型在农药计算方面的差异，包括：
1. 模型架构与算法差异
2. 农药运移过程处理方式
3. 输入参数要求
4. 输出结果格式与精度
5. 适用场景与局限性

## 模型背景资料

### SWAT+ 模型

**开发机构**: USDA Agricultural Research Service (USDA-ARS) + Texas A&M AgriLife Research

**模型类型**: 小流域到河流流域尺度的水质水量模拟模型

**主要功能**:
- 模拟地表水和地下水的质量与数量
- 预测土地利用、土地管理措施和气候变化的环境影响
- 广泛应用于土壤侵蚀控制、非点源污染控制和流域管理

**农药模块特点**:
- 支持农药降解、运移过程模拟
- 可模拟除草剂、杀虫剂、杀菌剂等多种农药类型
- 基于HRU（水文响应单元）进行分布式模拟
- 2025年最新技术笔记扩展了农药模拟能力

**应用案例**: 超过3000篇发表研究，约50篇专门针对农药模拟

**参考资料**:
- [SWAT+官方文档](https://swatplus.gitbook.io/docs)
- [SWAT+源代码文档](https://swat-model.github.io/swatplus)
- [Texas A&M SWAT](https://swat.tamu.edu)
- [Technical note: Extending SWAT2012 and SWAT+](https://hess.copernicus.org/articles/29/6703/2025/)
- [Pesticide transport simulation in tropical catchment](https://www.sciencedirect.com/science/article/abs/pii/S0269749114001511)
- [Advancing Pesticide Exposure Assessments with SWAT+](https://swat.tamu.edu/media/t4ijyo52/a3-1-swatplus_strassbourg_03jul2024.pdf)

### PWC 模型

**开发机构**: US EPA (Environmental Protection Agency)

**模型类型**: 农药水环境浓度计算工具

**主要功能**:
- 估算农药施用于土地表面后在水体中的环境浓度
- 用于农药风险评估和管理决策
- 支持地表水和地下水评估

**核心组件**:
- **PRZM3** (Pesticide Root Zone Model) - 模拟农药在根区的运移
- **VVWM** (Variable Volume Water Model) - 模拟水体中的农药浓度变化

**农药模块特点**:
- 基于场景的确定性模型
- 集成挥发、降解、吸附等过程
- 专注于监管风险评估应用
- 与PRZM3耦合使用

**应用场景**:
- EPA农药监管风险评估
- 濒危物种生物学评估
- 流域尺度农药运移建模

**参考资料**:
- [PWC用户手册](https://19january2021snapshot.epa.gov/sites/static/files/2015-12/documents/pwc_user_manual_12-8-15.pdf)
- [EPA农药风险评估模型](https://www.epa.gov/pesticide-science-and-assessing-pesticide-risks/models-pesticide-risk-assessment)
- [PWC用于土壤浓度估算](https://academic.oup.com/ieam/article/21/4/943/8043221)
- [PWC with PRZM](https://www.midatlanticrisa.org/data-tools/water-model-tool/items/pwc-using-the-przm.html)

## 项目结构

```
SWATvsPWC/
├── PWC_src/                    # PWC模型源代码
│   ├── Chemica_Transport.f90   # 化学物质运移
│   ├── Pesticide_Application.f90 # 农药施用
│   ├── Plant_Pesticide_Processes.f90 # 植物农药过程
│   ├── volatilization.f90      # 挥发过程
│   ├── TPEZ-WPEZ.f90          # 土壤-植物-水体运移
│   └── vvwm/                   # VVWM水体模块
│       ├── CoreCalculations.F90
│       ├── degradation.F90
│       └── VolumeAndWashout.F90
│
├── swatplus/                    # SWAT+模型源代码
│   ├── src/                    # Fortran源代码
│   ├── doc/                    # 文档
│   ├── refdata/                # 参考数据
│   └── test/                   # 测试用例
│
├── claude.md                    # 本规划文档
├── comparison/                  # 对比分析结果
│   ├── docs/                   # 对比文档
│   │   ├── swat_vs_pwc_presentation.md        # **Marp演示文稿（核心产出）**
│   │   ├── 01_model_comparison_report.md      # 综合对比报告
│   │   ├── 03_parameter_comparison_table.md   # 参数对照表
│   │   └── 04_application_guide.md             # 应用指南
│   └── analysis/               # 代码分析
│       └── 02_code_analysis_and_flowcharts.md  # 代码分析与流程图
└── docs/                        # 分析文档
```

## 研究计划

### 阶段一：模型代码分析

**SWAT+ 农药模块分析**:
1. 定位SWAT+源代码中的农药相关模块
2. 分析农药运移算法实现
3. 理解输入参数和输出变量
4. 识别关键计算公式

**PWC 农药模块分析**:
1. 分析`Pesticide_Application.f90`和`Chemical_Transport.f90`
2. 研究VVWM水体模块
3. 理解PRZM集成方式
4. 识别关键计算公式

### 阶段二：算法对比

对比以下方面：
1. **农药施用模拟**
   - 施用方式（叶面、土壤）
   - 施用时间和频率处理
   - 初始条件设定

2. **运移过程模拟**
   - 吸附/解吸过程
   - 降解过程（好氧/厌氧）
   - 挥发过程
   - 淋溶过程
   - 地表径流携带

3. **水文过程耦合**
   - 与水文模块的集成方式
   - 时间步长处理
   - 空间尺度处理

4. **输出结果**
   - 输出变量类型
   - 时间分辨率
   - 空间分辨率

### 阶段三：差异总结与应用建议

1. 总结两个模型的核心差异
2. 分析各自的优缺点
3. 提出模型选择建议
4. 识别可能的改进方向

## 预期成果

### 最终产出

**演示文稿**: [SWAT+ vs PWC 对比演示文稿](comparison/docs/swat_vs_pwc_presentation.md) - **项目核心目标**
- 基于源代码分析的模型对比PPT
- 涵盖模型概述、代谢产物模拟、物理过程、施用方法、挥发/吸附/数值方法、输出对比、场景对比等
- 使用 Marp 格式（基于 Warp 模板构建），支持导出为 PDF

### 已完成的文档

1. **技术文档**: [详细的模型对比分析报告](comparison/docs/01_model_comparison_report.md)
   - 模型概述与架构对比
   - 核心算法对比（水文耦合、数值方法）
   - 农药施用、降解、挥发、吸附过程对比
   - 输出结果对比
   - 模型优缺点分析

2. **代码分析**: [关键算法的代码注释和流程图](comparison/analysis/02_code_analysis_and_flowcharts.md)
   - SWAT+模块架构分析
   - PWC模块架构分析
   - 算法流程图对比
   - 关键函数详解

3. **参数对比**: [输入参数要求的对照表](comparison/docs/03_parameter_comparison_table.md)
   - 农药物理化学参数
   - 施用参数
   - 降解参数
   - 吸附参数
   - 挥发参数
   - 植物相关参数
   - 水文与环境参数
   - 数值计算参数

4. **应用指南**: [针对不同场景的模型选择建议](comparison/docs/04_application_guide.md)
   - 快速决策指南
   - 场景选择矩阵
   - SWAT+使用指南
   - PWC使用指南
   - 数据准备清单
   - 模型校准建议
   - 常见问题解答

## 参考文献来源

### SWAT+ 相关
- [Extending SWAT2012 and SWAT+ models (HESS, 2025)](https://hess.copernicus.org/articles/29/6703/2025/)
- [Advancing Pesticide Exposure Assessments with SWAT+](https://swat.tamu.edu/media/t4ijyo52/a3-1-swatplus_strassbourg_03jul2024.pdf)
- [SWAT+ Watershed Simulation of Wetlands and Pesticide](https://ag.purdue.edu/department/arge/industry/wetlands/_media-wetlands/arnold_wetland_conf_purdue2020-autosaved-a.pdf)
- [Simulation of Pesticide and Metabolite Concentrations](https://www.mdpi.com/2073-4441/14/9/1332)
- [A Review of Pesticide Fate and Transport Simulation](https://pubmed.ncbi.nlm.nih.gov/30884273/)

### PWC 相关
- [PWC User Manual (EPA, 2015)](https://19january2021snapshot.epa.gov/sites/static/files/2015-12/documents/pwc_user_manual_12-8-15.pdf)
- [Models for Pesticide Risk Assessment (EPA)](https://www.epa.gov/pesticide-science-and-assessing-pesticide-risks/models-pesticide-risk-assessment)
- [PWC tool for soil concentration estimates](https://academic.oup.com/ieam/article/21/4/943/8043221)
- [PWC using PRZM](https://www.midatlanticrisa.org/data-tools/water-model-tool/items/pwc-using-the-przm.html)

---

*文档创建时间: 2026-02-01*
