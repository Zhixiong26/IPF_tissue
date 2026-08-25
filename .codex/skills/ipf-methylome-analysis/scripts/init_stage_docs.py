#!/usr/bin/env python3
"""Create non-destructive README.md and Report.md templates for stage scripts."""
import argparse
from pathlib import Path


README_TEMPLATE = """# {stage}

本目录保存 {stage} 阶段的版本化分析脚本和调度入口。

This directory contains versioned analysis scripts and scheduler entry points for the {stage} stage.

## 目录约定 / Directory contract

- 执行前记录所需输入和元数据。 / Document required inputs and metadata before execution.
- 分析产物写入 `Results/{stage}/`。 / Write analysis products below `Results/{stage}/`.
- 运行及调度日志写入 `Scripts/{stage}/logs/`。 / Write runtime and scheduler output below `Scripts/{stage}/logs/`.
- 添加脚本时同步记录环境、命令、参数和断点续跑行为。 / Record environments, commands, parameters, and restart behavior as scripts are added.

## 执行入口 / Entry points

尚未登记可执行流程。 / No executable workflow has been registered yet.
"""


REPORT_TEMPLATE = """# {stage} 分析报告 / Analysis report

状态 / Status: 未开始 / not started

## 目标 / Objective

记录本阶段的具体问题和完成标准。 / Record the concrete question and completion criteria for this stage.

## 输入与配置 / Inputs and configuration

记录输入清单、软件/环境版本、参数和调度资源。 / Record input manifests, software/environment versions, parameters, and scheduler resources.

## 结果与质控 / Results and QC

每次运行后记录已验证的输出和定量质控；任务提交不等于结果完成。 / Record verified outputs and quantitative QC after each run; a submitted job is not a completed result.

## 决策与限制 / Decisions and limitations

记录未决选项、失败尝试、流程偏差和解释限制。 / Record unresolved choices, failed attempts, deviations, and interpretation limits.
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("project_dir", type=Path)
    parser.add_argument("stages", nargs="+")
    args = parser.parse_args()

    for stage in args.stages:
        stage_dir = args.project_dir / "Scripts" / stage
        stage_dir.mkdir(parents=True, exist_ok=True)
        for filename, template in (
            ("README.md", README_TEMPLATE),
            ("Report.md", REPORT_TEMPLATE),
        ):
            path = stage_dir / filename
            if path.exists():
                print("preserved {}".format(path))
                continue
            path.write_text(template.format(stage=stage), encoding="utf-8")
            print("created {}".format(path))


if __name__ == "__main__":
    main()
