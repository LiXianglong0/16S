# 方案 B：日常推 Gitee，自动同步到 GitHub

目标：本地只推 **Gitee**（国内稳）；**GitHub** 由流水线自动跟上。

- Gitee（日常，已建好）：`https://gitee.com/git_lixianglong/16S`
- GitHub（镜像）：`https://github.com/LiXianglong0/16S`

---

## 一、一次性配置（约 10 分钟）

### 1. 在 Gitee 新建空仓库

1. 打开 https://gitee.com → 新建仓库  
2. 仓库名：`16S`  
3. **不要**勾选「使用 Readme 文件初始化」  
4. 创建后复制 HTTPS 地址，形如：  
   `https://gitee.com/你的用户名/16S.git`

### 2. 把本地 `origin` 改成 Gitee（本机执行）

把下面的 `你的用户名` 换成真实 Gitee 用户名：

```bash
cd /mnt/d/study/16S   # 或 Windows：cd D:\study\16S

# 原 GitHub 留作备用远程名 github
git remote rename origin github

# 日常推送用 origin → Gitee
git remote add origin https://gitee.com/你的用户名/16S.git

git remote -v
# 应看到：
# origin  → gitee
# github  → github.com/LiXianglong0/16S.git
```

也可直接运行（先改脚本里的用户名）：

```bash
bash scripts/setup-gitee-origin.sh 你的Gitee用户名
```

### 3. 首次推到 Gitee

```bash
git push -u origin main
```

### 4. 配置「Gitee → GitHub」自动同步（二选一）

#### 方式 A（推荐）：GitHub Actions 定时从 Gitee 拉取

在 **GitHub** 仓库添加 Secrets（Settings → Secrets and variables → Actions）：

| Name | 值 |
|---|---|
| `GITEE_REPO` | `https://gitee.com/你的用户名/16S.git` |
| `GITEE_TOKEN` | Gitee 私人令牌（仓库若公开可留空，私有必填） |

Gitee 令牌：Gitee → 设置 → 安全设置 → 私人令牌 → 勾选 `projects`。

本仓库已包含工作流文件：

`.github/workflows/sync-from-gitee.yml`

- 每小时自动跑一次，也可在 Actions 页手动 **Run workflow**  
- **前提**：该 yml 至少成功出现在 GitHub 仓库里一次  

若本地仍推不上 GitHub，用下面「方式 B」先把代码放到 GitHub，或在 GitHub 网页：  
**Creating new file** → 路径填 `.github/workflows/sync-from-gitee.yml` → 把本仓库同名文件内容粘贴上去 → Commit。

#### 方式 B：Gitee 网页镜像（有则优先）

Gitee 仓库 → **管理** → **仓库镜像设置 / 同步仓库**（名称因套餐略有不同）：

- 镜像地址：`https://github.com/LiXianglong0/16S.git`  
- 用户名：GitHub 用户名  
- 密码：GitHub **Personal Access Token**（勾选 `repo`）  
- 开启推送同步 / 定时同步  

---

## 二、以后日常怎么用

```bash
git add .
git commit -m "说明本次改动"
git push origin main          # 只推 Gitee
```

等 Actions / 镜像跑完后，GitHub 的 `main` 会自动更新。  
偶尔想手动推 GitHub（网络好时）：

```bash
git push github main
```

---

## 三、检查是否同步成功

1. Gitee 仓库提交是否最新  
2. GitHub → Actions → **Sync from Gitee** 是否绿勾  
3. GitHub 仓库代码是否与 Gitee 一致  

---

## 四、注意

- Issues / PR 不会互通，只同步提交历史  
- `bak/`、数据库大文件、原始 FASTQ 已在 `.gitignore` 中忽略，网盘另行存放  
- 不要把 Token 写进代码，只放在 GitHub Secrets / Gitee 密钥里  
