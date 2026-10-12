# 从 omicsapp:1.0 升级：日志、签名密钥、端口、备份 —— 操作清单

适用于还在运行 `omicsapp:1.0`（按旧版部署文档手动打标签的镜像）的服务器。按顺序执行，每一步验证通过再做下一步。

| 本清单做的事 | 对应 |
|---|---|
| 每晚带日期的快照备份 + 每周自动恢复演练 | 第 1 步 |
| 每个用户会话的 R 日志落盘；ShinyProxy 自身日志自动轮转 | 第 2 步 |
| 项目文件签名密钥；容器内不再尝试联网更新基因集 | 第 3 步 |
| 用户会话端口只对本机开放 | 第 4 步 |

- **预计时间：** 第 1 步约 30 分钟（取决于数据量）；第 2–6 步约 20 分钟。
- **停机：** 只有第 5 步重启 ShinyProxy 一次，在线会话会断开；应用有自动保存，用户刷新后可恢复。如果同时要升级到新镜像，把两件事放在同一次重启里做（见第 5 步）。
- **独立于镜像版本：** 第 1、2、4 步对 1.0 镜像同样有效。第 3 步的两项环境变量旧镜像会忽略，换成新镜像后才生效。
- **详细原理：** 见 [`README.md`](README.md) 的 "Operations"、"Project files are signed"、"Where the logs are" 各节。

以下命令中：

- `~/omicsApp` 是你从工作电脑 rsync 过来的代码（[`README.md`](README.md) 第 0 步）；
- `/opt/omicsApp` 是供 root 定时任务使用的部署副本（第 0.2 步创建）；
- `$EDITOR` 换成你习惯的编辑器（如 `nano`、`vim`）。

---

## 第 0 步：准备（不影响服务）

### 0.1 同步代码

在工作电脑上：

```bash
cd ~/SciProject/CHISSS/omicsApp
git fetch origin && git checkout <要部署的分支> && git pull
git status --short                 # 必须为空
git rev-parse HEAD > REVISION
rsync -avz --delete --exclude-from=deploy/rsync.exclude ./ <用户名>@<服务器>:~/omicsApp/
```

### 0.2 建立 root 所有的部署副本

```bash
sudo rsync -a --delete --chown=root:root ~/omicsApp/ /opt/omicsApp/
```

备份和恢复演练由 root 的 cron 运行，脚本位置是 `/opt/omicsApp/deploy/scripts/`。不要让 root 直接运行你家目录里的脚本：任何能改写 `~/omicsApp` 的人都能借此以 root 身份执行命令。

以后每次同步代码后，都要重复执行这一行。

### 0.3 备份当前的 ShinyProxy 配置

```bash
sudo cp -p /etc/shinyproxy/application.yml /etc/shinyproxy/application.yml.bak-$(date +%F)
```

### 0.4 确认备份盘和数据盘不是同一块物理盘

```bash
df -h /srv /backup     # 两行的 Filesystem 列应当是不同的设备
```

如果 `/backup` 不存在，第 1.1 步会创建它。创建后再执行一次这条命令确认：它应该落在系统盘（SSD），而用户数据在 `/srv`（HDD）。

**如果两者是同一块盘，先停下来。** 这种情况下备份和数据会同时损坏，需要另外找一块盘或另一台机器。

### 0.5 生成签名密钥，并单独保存一份

```bash
openssl rand -hex 32
```

把输出的 64 位十六进制字符串保存到密码管理器（或其他不在这台服务器上的安全位置），第 3 步要用。

- **这个密钥丢了，所有用户的项目文件都打不开。**
- 密钥会随第 1 步的备份一起保存，但只在服务器上有一份是不够的。

---

## 第 1 步：备份与恢复演练

先做备份，后面每一步改动都有退路。

### 1.1 配置

```bash
sudo mkdir -p /etc/omicsapp /backup/omicsapp
sudo cp /opt/omicsApp/deploy/backup.env.template /etc/omicsapp/backup.env
sudo chmod 600 /etc/omicsapp/backup.env
sudo $EDITOR /etc/omicsapp/backup.env
```

至少确认这几行：

| 设置 | 填什么 |
|---|---|
| `BACKUP_ROOT` | `/backup/omicsapp`（与数据不同的物理盘） |
| `REPO_DIR` | `/opt/omicsApp` |
| `BACKUP_REMOTE` | 另一台机器，如 `backup@另一台主机:/srv/backup/omicsapp`。留空则只有本机一份，机器损坏、被盗或重装系统会同时失去数据和备份 |
| `ALERT_WEBHOOK` / `ALERT_EMAIL` | 备份失败时通知谁（钉钉、企业微信、Slack 等接受 `{"text": "..."}` 的 webhook，或邮箱） |
| `KEEP_DAILY` / `KEEP_WEEKLY` | 默认保留 14 天的每日快照和 8 周的每周快照。空间紧张时可调小 |

### 1.2 设置了 `BACKUP_REMOTE` 时：配置 root 的 ssh 免密登录

cron 无法回答密码提示，所以必须用密钥登录：

```bash
sudo ssh-keygen -t ed25519 -f /root/.ssh/omicsapp_backup -N ''
sudo ssh-copy-id -i /root/.ssh/omicsapp_backup.pub <远程用户>@<远程主机>
```

然后在 `/root/.ssh/config` 里为该主机指定这个密钥：

```
Host <远程主机>
    IdentityFile /root/.ssh/omicsapp_backup
```

最后确认免密登录可用：

```bash
sudo ssh <远程用户>@<远程主机> true    # 不应提示输入密码
```

### 1.3 手动运行一次

```bash
sudo /opt/omicsApp/deploy/scripts/backup.sh
```

- **预期：** 最后一行以 `done` 结尾。
- **第一次运行：** 会复制全部用户数据，耗时取决于数据量。之后每晚只复制有变化的文件。

检查结果：

```bash
sudo ls /backup/omicsapp/snapshots/                    # 有一个带日期的目录
sudo du -sh /backup/omicsapp                           # 实际占用（硬链接只计一次）
sudo ls /backup/omicsapp/latest/config/host/etc/shinyproxy/   # 有 application.yml
```

最后一行确认线上的 ShinyProxy 配置也被备份了：第 3 步加的签名密钥就保存在这个文件里。

### 1.4 安装定时任务

```bash
sudo cat /etc/cron.d/omicsapp-backup 2>/dev/null   # 旧版本的定时任务（如果有），下一行会替换它
sudo cp /opt/omicsApp/deploy/cron/omicsapp-backup /etc/cron.d/omicsapp-backup
sudo chmod 644 /etc/cron.d/omicsapp-backup
sudo cp /opt/omicsApp/deploy/logrotate/omicsapp /etc/logrotate.d/omicsapp
```

- **定时：** 每晚 02:20 备份，每周日 04:10 恢复演练，输出写到 `/var/log/omicsapp-backup.log`。
- **为什么必须 `chmod 644`：** cron 会静默忽略组或其他用户可写的文件。权限不对时不会报错，也不会有备份。

### 1.5 手动运行一次恢复演练

```bash
sudo /opt/omicsApp/deploy/scripts/restore_check.sh
```

**预期：** 以 `all checks passed` 结尾。它会做这些检查：

- 最新备份是否在 36 小时以内；
- 校验值是否一致；
- 把账号数据库恢复到一个临时 Postgres 容器，并核对账号数；
- 随机打开一个备份里的项目；
- 设置了远程备份时，远程是否也有最新快照。

不影响正在运行的服务。

### 1.6 旧备份

旧版本的定时任务留下的 `/backup/omicsapp/keycloak-*.sql.gz` 等文件，可以在新的快照连续正常运行一周后删除，以释放空间：

```bash
sudo ls -la /backup/omicsapp/      # 只删旧的 keycloak-?.sql.gz；不要动 snapshots/ 和 latest
```

---

## 第 2–4 步：修改 ShinyProxy 配置

三步都改 `/etc/shinyproxy/application.yml`，最后在第 5 步重启一次。

**直接编辑服务器上这份文件，不要用 `deploy/shinyproxy/application.yml` 覆盖它。** 线上文件里有 OIDC client secret，模板里没有。

### 2.1 先创建日志目录

**必须在改配置之前创建**，否则 ShinyProxy 启动时找不到可写的目录：

```bash
sudo install -d -o shinyproxy -g shinyproxy -m 750 /var/log/shinyproxy/containers
```

（logrotate 的轮转规则已在 1.4 步安装。）

### 2.2–4 编辑配置

```bash
sudo $EDITOR /etc/shinyproxy/application.yml
```

按下面的位置加入标 `# 新增` 的行。YAML 对缩进敏感，**用空格，不要用 Tab**，并与周围已有的行对齐。

```yaml
proxy:
  # ...（原有内容不变）...
  default-stop-proxy-on-logout: true

  container-log-path: /var/log/shinyproxy/containers      # 新增（第 2 步）：每个会话的 R 日志

  docker:
    internal-networking: false
    target-bind-ip: 127.0.0.1                               # 新增（第 4 步）：会话端口只对本机开放

  specs:
    - id: omicsapp
      # ...（原有内容不变）...
      container-env:
        OMICSAPP_QUOTA_GB: "50"
        OMICSCORE_GENESET_TTL_DAYS: "0"                     # 新增（第 3 步）：容器内不联网更新基因集
        OMICSAPP_SIGNING_KEY: <第 0.5 步生成的密钥>          # 新增（第 3 步）：项目文件签名密钥

logging:
  file:
    name: /var/log/shinyproxy/shinyproxy.log
  logback:                                                  # 新增（第 2 步）：ShinyProxy 自身日志轮转
    rollingpolicy:
      max-file-size: 10MB
      max-history: 30
      total-size-cap: 1GB
```

保存后确认文件属主（ShinyProxy 以 `shinyproxy` 用户运行，读不到就会启动失败）：

```bash
sudo chown shinyproxy:shinyproxy /etc/shinyproxy/application.yml
sudo chmod 600 /etc/shinyproxy/application.yml
```

关于签名密钥（第 3 步）：

- **以后不要删除这一行。** 删除后每个用户的目录会各自生成新密钥，然后拒绝打开用这个密钥签名的全部项目。要恢复访问，只能把原密钥重新加回来。
- **第一次被新镜像使用时：** 每个用户目录会对已有项目做一次结构检查并签名，项目内容本身不会被改写。未通过检查的文件会被拒绝打开，并写进 `/var/log/shinyproxy/containers/` 下的日志。
- **配合第 1 步的备份：** 这一行改完后，下一次备份会把含密钥的配置一起保存。也可以立即手动运行一次 `sudo /opt/omicsApp/deploy/scripts/backup.sh`。

---

## 第 5 步：重启 ShinyProxy（唯一一次停机）

提前通知在线用户，然后二选一。

**同时升级到新镜像**（已按 [`README.md`](README.md) 用 `build_image.sh` 构建并试运行过）：

```bash
sudo /opt/omicsApp/deploy/scripts/rollback.sh omicsapp:0.2.0-<提交号>
```

它会改写镜像设置、重启 ShinyProxy，并等待它恢复响应。

**只应用本清单、暂不换镜像：**

```bash
sudo systemctl restart shinyproxy
```

确认已正常启动：

```bash
sudo systemctl status shinyproxy --no-pager       # active (running)
sudo tail -n 50 /var/log/shinyproxy/shinyproxy.log # 没有 ERROR，最后有 "Started ..."
```

如果起不来，最常见的原因是 YAML 缩进写错或文件属主不对。把配置恢复成 0.3 步的备份再重启，然后检查改动：

```bash
sudo cp -p /etc/shinyproxy/application.yml.bak-<日期> /etc/shinyproxy/application.yml
sudo systemctl restart shinyproxy
```

---

## 第 6 步：验证

登录应用，打开一个已有项目，再新建或保存一次，然后逐项检查：

| 检查 | 命令 | 预期 |
|---|---|---|
| 会话日志落盘 | `sudo ls -la /var/log/shinyproxy/containers/` | 刚才的会话有一个 `.log` 文件 |
| 端口只对本机 | 会话进行中运行 `sudo ss -tlnp \| grep ':2[0-9]\{4\} '` | 地址都是 `127.0.0.1`，没有 `0.0.0.0` 或 `*` |
| 签名已启用（新镜像） | `sudo ls -a /srv/omicsapp/users/<用户目录>/` | 有 `.omicsapp-signatures`；**没有** `.omicsapp-signing-key`（用的是统一密钥） |
| 已有项目都能打开 | `sudo grep -l "was not saved by this app\|has been changed since" /var/log/shinyproxy/containers/*.log` | 没有输出 |
| 备份含线上配置 | `sudo grep -c OMICSAPP_SIGNING_KEY /backup/omicsapp/latest/config/host/etc/shinyproxy/application.yml`（第 3 步之后的备份） | `1` |
| ShinyProxy 日志轮转 | 一天后：`sudo ls /var/log/shinyproxy/` | 出现按日期归档的旧日志 |
| 定时备份在跑 | 第二天：`sudo tail /var/log/omicsapp-backup.log` | 有当晚的运行记录，结尾为 `done` |

用户目录名是 UUID，要知道哪个目录属于谁，用：

```bash
sudo /opt/omicsApp/deploy/scripts/list_users.sh
```

---

## 回退

| 要撤销的 | 做法 |
|---|---|
| 新镜像 | `sudo /opt/omicsApp/deploy/scripts/rollback.sh omicsapp:1.0` |
| 第 2、4 步的配置 | 从 `application.yml` 删除对应行，然后 `sudo systemctl restart shinyproxy` |
| 第 3 步的签名密钥 | **用户用新镜像打开过项目之后，不要删除。** 删除只会让所有项目打不开 |
| 第 1 步的定时备份 | `sudo rm /etc/cron.d/omicsapp-backup`（已有快照保留在 `/backup/omicsapp/`） |

回退到 1.0 镜像后又再次升级时，先执行一次：

```bash
sudo rm /srv/omicsapp/users/*/.omicsapp-signatures
```

这样新镜像会重新检查并签名用户在旧版本期间保存的项目。
