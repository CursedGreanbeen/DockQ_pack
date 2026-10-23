# DockQ — пакетная оценка метрик DockQ для моделей антитело-антигенных комплексов, полученных после фолдинга с AlphaFold

Считает DockQ для каждой предсказанной модели (эталон - экспериментальная структура): одна модель —
один вызов CLI DockQ как отдельный процесс (`run.py`), разбор вывода — в `analyze.py`._interfaces_
помечаются как AB–AG по ANARCI-аннотациям цепей (необходимо предварительно проаннотировать!); результаты построчно дописываются в TSV-отчёты.
Рассчитан как на стандартный фолдинг, так и на multiseed-раскладку (10 комплексов × 500 сидов × 5 сэмплов) с опциональной
предфильтрацией по ranking CSV.

## Структура и конвейер

Конвейер (единственная точка входа — `cli.py`):
**ranking CSV (опц. предфильтр) → поиск моделей → поиск эталонов → запуск DockQ → парсинг stdout → TSV-отчёты**

| Файл | Роль |
|---|---|
| `cli.py` | Основной скрипт: поиск моделей, предфильтр, параллельный запуск, запись отчётов |
| `ranking.py` | `filter_models_by_ranking()` — топ-N сэмплов на сид из `<pdb>/<pdb>_ranking_scores.csv` |
| `run.py` | Обёртка над CLI DockQ (бинарник захардкожен: `/data/user/mullagaliamova/.conda/envs/dockq_env/bin/DockQ`) |
| `analyze.py` | Парсинг stdout DockQ; проверка AB–AG; бины качества (High ≥0.80, Medium ≥0.49, Acceptable ≥0.23, иначе Incorrect); агрегация |
| `mapping.py` | Загрузка ANARCI TSV → наборы AB/AG-цепей на PDB (`ANTIBODY`→AB, `NOT_RECOGNIZED`/`NO_ANTIBODY_DOMAIN`→AG) |
| `resolution.py` | Пост-обработка: добавление колонки resolution из SAbDab к TSV DockQ |
| `recreate_interface_report.py` | Пересборка AB–AG-отчёта из detailed-отчёта (цепи в namespace модели) |
| `recreate_interface_report_ref.py` | То же, но фильтр в **reference**-namespace через `anarci_report_ref.tsv` (корректный вариант) |
| `__init__.py` | Экспорты пакета |

Структура данных (multiseed): модели — `<task-dir>/<pdb>/seed-<seed>_sample-<k>/*_model.cif`,
ranking_score — `<task-dir>/<pdb>/<pdb>_ranking_scores.csv` (колонки `seed,sample,ranking_score`),
Эталоны — `blasted-cif/<pdb>.cif`. `find_model_files()` берёт по одной модели из каждого seed-каталога
и пропускает верхнеуровневый `<pdb>_model.cif`; для комплексов без seed-каталогов возвращается
единственная верхнеуровневая модель.

## Запуск

Работать из каталога scripts как cwd; интерпретатор — `/opt/conda/bin/python3.10`.
`cli.py`/`recreate_*` сами добавляют родительский каталог пакета в `sys.path`, поэтому абсолютные
пути к скриптам работают:

```bash
cd /mnt/1857e392-c689-44a8-9c3e-e7f9beab82f6/scripts
/opt/conda/bin/python3.10 dockq/cli.py ...
```

### cli.py (основной батч)

```bash
/opt/conda/bin/python3.10 dockq/cli.py \
  --s3data-dir /mnt/f0771a75-815c-4553-ba74-9aa70f57f995/<task-dir> \
  --cifs-dir /mnt/1857e392-c689-44a8-9c3e-e7f9beab82f6/blasted-cif \
  --top-per-seed 1 \
  --interface-report /mnt/1857e392-c689-44a8-9c3e-e7f9beab82f6/results/reports/dockq_500seeds_interfaces.tsv \
  -o /mnt/1857e392-c689-44a8-9c3e-e7f9beab82f6/results/reports/dockq_500seeds_detailed.tsv
```

Ключевые флаги: `--s3data-dir`, `--cifs-dir`, `--cifs-filtered-dir`
(поиск эталона: `{pdb}_cropped.cif`, `{pdb}.cif` — сначала в filtered, затем в all-каталоге),
`-o/--output` (detailed TSV по всем интерфейсам), `--interface-report` (только AB–AG),
`--top-per-seed N` (предфильтр по ranking_score — **подключён**: каждый CSV читается один раз и кэшируется;
при отсутствии CSV/сида модель сохраняется с предупреждением), `--pdb-id`, `--limit N`,
`--workers N` (по умолчанию 4), `--log-dir`, `--anarci-report` (по умолчанию
`results/reports/anarci_report.tsv`), `--no-chain-mapping`, `--summary`,
`--group-by pdb_id|ab_chains|ab_ag_pair`.

Подводные камни: built-in-дефолты (`--s3data-dir`, `--cifs-*`) указывают на легаси-пути
`/home/mullagaliamova/...` — флаги каталогов передавать обязательно. Флаг `--summary` принимается,
но **не работает**: `main()` никогда не вызывает `write_summary_report()`, файл сводки не пишется.
DockQ запускается без `--mapping`; AB/AG-фильтрация — постфактум по `native_chains` (см.
`chain_mapping_findings.md`: namespace'ы цепей эталона и модели различаются — использовать
`recreate_interface_report_ref.py`).

### resolution.py

```bash
/opt/conda/bin/python3.10 dockq/resolution.py \
  --sabdab /path/to/sabdab_2024_hr_full_final.tsv \
  --dockq /mnt/1857e392-c689-44a8-9c3e-e7f9beab82f6/results/reports/dockq_500seeds_detailed.tsv \
  -o /mnt/1857e392-c689-44a8-9c3e-e7f9beab82f6/results/reports/dockq_detailed_with_resolution.tsv
```

Входы: `--sabdab` (по умолчанию — легаси-TSV в `/home/...`), `--dockq` (обязателен, нужна колонка
`pdb_id`). Выход: входной TSV + колонка `resolution` (путь по умолчанию — легаси-каталог отчётов).

### recreate_interface_report.py / recreate_interface_report_ref.py

Позиционные аргументы, без argparse:

```bash
/opt/conda/bin/python3.10 dockq/recreate_interface_report_ref.py \
  <detailed.tsv> <interface.tsv> [anarci_ref.tsv]   # по умолчанию: results/reports/anarci_report_ref.tsv
```
На случай непредвиденных ошибок и остановки скрипта на середине работы:
`recreate_interface_report.py <detailed.tsv> <interface.tsv>` фильтрует по колонкам
`ab_chains`/`ag_chains` из detailed-отчёта (namespace модели); вариант `_ref` фильтрует по
`anarci_report_ref.tsv` и переписывает эти колонки в reference-namespace.
