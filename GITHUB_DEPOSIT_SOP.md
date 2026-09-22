# GitHub / Zenodo 复现包存放操作手册（中文）

> 适用：每个科研项目配一个复现包仓库（kebab-case 主题名）。本手册为 `pod-epigenetic-immune-repositioning` 的存放流程。

## 一、GitHub 存放

1. 本地仓库已建好（`pod-epigenetic-immune-repositioning/`），含 README / LICENSE / CITATION.cff / 分析脚本。
2. 初始化并推送：
   ```bash
   cd pod-epigenetic-immune-repositioning
   git init
   git add .
   git commit -m "init: reproducibility package for POD epigenetic immune repositioning"
   gh repo create yyx-4113/pod-epigenetic-immune-repositioning --public --source=. --push
   ```
3. 打版本标签（与稿件版本一致，如 v1.0.0）：
   ```bash
   git tag v1.0.0 && git push origin v1.0.0
   ```
   GitHub Actions（`.github/workflows/release.yml`）会自动生成 Release 并打包分析脚本。

## 二、Zenodo 归档（获取 DOI，用于数据可用性声明）

1. 用 GitHub 账号登录 Zenodo → "New upload" → "GitHub" → 选择本仓库与对应 tag。
2. 填写：标题、作者（Yongxin Yang, ORCID 0009-0004-9698-6552）、许可证 MIT、关键词。
3. 发布后复制 Zenodo DOI，写入稿件 Data Availability Statement（**实名仓库 URL + Zenodo DOI，禁止 "available on request"**）。

## 三、稿件中引用格式（Vancouver）

> Yang Y. Postoperative delirium epigenetic immune hub and in-silico drug repositioning — reproducibility package [Internet]. GitHub; 2026. Available from: https://github.com/yyx-4113/pod-epigenetic-immune-repositioning. DOI: 10.5281/zenodo.XXXXXXX.

## 四、更新纪律

- 每次分析有重大改动 → 提交并打新 tag → 同步 Zenodo 新版本。
- 脚本内软件版本（R/Bioconductor/Python）在投稿前核对并钉死到 README「Software versions」。
- 作者头衔只写姓名，不填 MD/PhD（无研究生学位）。
