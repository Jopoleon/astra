# 01. ASTRA 8: обзор кода, история и карта исходников

**О чём глава.** Что такое ASTRA, откуда она взялась и чем ASTRA 8 отличается от ASTRA 7; кто её поддерживает и под какой лицензией; полная карта дерева исходников `astra-src` (каждый каталог и файл верхнего уровня); идея «репозиторная область / пользовательская область»; общая схема потока данных от модели и эксперимента до результата; что концептуально вычисляет один расчёт. Детали сборки и запуска — в [02-build-and-run.md](02-build-and-run.md), остальное — в профильных главах (ссылки по тексту).

Состояние: исходники `astra-src`, MPCDF `67281112` (2026-10-01; совпадает с `upstream/main`). Рабочая установка для сверки — `~/astra/a8` (коммит `e66f6ba1`, отличается от `67281112` двумя строками: путь SPIDER в `exe/astra_rc` и `src/misc/imas_ids.f90`).

Контекст монорепозитория (установщик, проблемы Ubuntu, форк) здесь не повторяется:

- [`../../AGENTS.md`](../../AGENTS.md) — сводка известного об ASTRA 8 и правила работы;
- [`../../installer/docs/DETAILED.md`](../../installer/docs/DETAILED.md) — как устроена сборка, проблемы на Ubuntu и их решения, замеры;
- [`../../docs/astra/platforms.md`](../../docs/astra/platforms.md) — площадки `platform/env.*`;
- [`../../docs/local-fork.md`](../../docs/local-fork.md) — устройство git в `astra-src/` и журнал правок.

## Содержание

1. [Что такое ASTRA](#1-что-такое-astra)
2. [История: от ASTRA 1 до ASTRA 8.6](#2-история-от-astra-1-до-astra-86)
3. [Что изменилось в ASTRA 8 по сравнению с ASTRA 7](#3-что-изменилось-в-astra-8-по-сравнению-с-astra-7)
4. [Авторы, сопровождающие, лицензия](#4-авторы-сопровождающие-лицензия)
5. [Репозиторная и пользовательская области](#5-репозиторная-и-пользовательская-области)
6. [Карта дерева исходников](#6-карта-дерева-исходников)
7. [Документация в `docu/`](#7-документация-в-docu)
8. [Поток данных: от модели до CDF](#8-поток-данных-от-модели-до-cdf)
9. [Что вычисляет один расчёт](#9-что-вычисляет-один-расчёт)
10. [Расхождения документации и кода](#10-расхождения-документации-и-кода)
11. [Где смотреть в исходниках](#11-где-смотреть-в-исходниках)

---

## 1. Что такое ASTRA

**ASTRA** — Automated System for TRansport Analysis («автоматизированная система транспортного анализа»). Это 1.5-мерный транспортный код для токамаков (и, с ограничениями, стеллараторов): по радиусу решаются одномерные уравнения переноса частиц, энергии, тока (полоидального потока) и вращения, усреднённые по магнитным поверхностям, а форма поверхностей берётся из двумерного равновесия. Отсюда «1.5D».

Ключевая особенность, которая отличает ASTRA от «обычной» программы и сохранилась с 1980-х: **это не готовое приложение, а генератор программ**. Цитата из основного руководства (`docu/section12.tex`, раздел «Overview of the Astra code»):

> First of all, Astra is not a ready-to-run application, but a tool for building customized computer codes for the solution of a variety of transport problems in magnetically confined plasmas.

Пользователь описывает транспортную задачу на собственном языке ASTRA (файл модели `equ/<имя>`): какие уравнения решать, какие коэффициенты переноса и источники брать, что рисовать. Система превращает описание в Fortran-код, компилирует его вместе с библиотеками ядра и запускает. В ASTRA 8 эту роль выполняет Python-парсер `pyparse/`, а сборку — `make` (подробно — [02-build-and-run.md](02-build-and-run.md) и [03-model-language.md](03-model-language.md)).

Состав системы по `docu/section12.tex` (раздел «Introduction»), с соответствием в коде ASTRA 8:

| компонент по руководству | где в ASTRA 8 |
|---|---|
| обширная библиотека модулей (физические процессы, обработка данных) | `src/astra`, `src/feqis`, `src/nbi`, `src/misc`; пользовательские `sbr/`, `fnc/`, `fml/` (из `repo_user_area/`) |
| «надзирающая оболочка», которая отслеживает изменения библиотек и собирает приложение | `exe/as_exe`, `exe/as_exe.py`, `exe/Build`, `exe/Makefile`, `pyparse/` |
| графический интерфейс управления расчётом и просмотра данных | `src/graph` (X11), см. [08-graphics-and-output.md](08-graphics-and-output.md) |
| интерфейс к экспериментальной базе данных | файлы эксперимента `exp/`, U-файлы `udb/`, `pyparse/exp_parser.py`, `pyparse/ufiles.py`, `src/astra/read_input.f90`; см. [04-experiment-data.md](04-experiment-data.md) |
| встроенная справка и примеры | примеры моделей `repo_user_area/equ/`, данные `repo_user_area/exp/`, `udb/`; встроенной справки в A8 не найдено |

Два режима работы, описанные ещё в руководстве 2002 года, сохранились:

- **активный (интерактивный)** — окно X11, пользователь видит эволюцию профилей, может останавливать расчёт и менять параметры на лету (`TASK = "RUN"`);
- **фоновый (batch)** — без окна, результаты сохраняются и смотрятся потом (`TASK = "BGD"`, флаг `-batch`).

## 2. История: от ASTRA 1 до ASTRA 8.6

Хронология собрана из `docu/section12.tex`, `docu/addendum_EF.tex`, `docu/summary_of_changes_a70.tex`, `docu/addendum_A8.tex`, `docu/IPC.tex`, `repo_user_area/exe/version` и истории git.

| период | версия | что известно | источник |
|---|---|---|---|
| конец 1960-х | — | первые транспортные коды (Кадомцев–Погуце; Днестровский–Костомаров). Г. В. Переверзев был среди пользователей кода Днестровского–Костомарова | `docu/section12.tex`, Introduction |
| конец 1980-х | первая ASTRA | автоматический построитель кода, Курчатовский институт (Москва) | `docu/section12.tex` |
| с 3.0 | 3.x | возможности одномерного моделирования стеллараторов | `docu/section12.tex` |
| февраль 2002 | 5.3 | основное руководство IPP 5/98 (Pereverzev, Yushmanov): Sun, IBM, DEC UNIX, IBM PC под Red Hat; C shell, Fortran/C, X11; «~20 МБ на диске» | `docu/main.tex`, `docu/section12.tex` |
| — | 6.1 | параллельные процессы через IPC (до 20), вызов модуля с `&` в модели | `docu/IPC.tex` |
| — | 6.2.1 | неявная равнораспределённость (equipartition) Te/Ti, IPC, схема для «жёсткого» (stiff) транспорта — в коде есть, в руководстве нет | `docu/summary_of_changes_a70.tex` |
| 2012–2013 | 7.0 | E. Fable, «Addendum with changes and new features»: нормированная сетка \(x=\sqrt{\Phi/\Phi_b}\), уравнение тороидального момента, устойчивая связь с равновесием (P′, FF′), интерфейсы TORBEAM, TGLF, QuaLiKiZ, TORIC, STRAHL, NEO, restart-файлы, блочные комментарии `%`, выбор блоков `#`, Ньютоновский решатель для стационара | `docu/addendum_EF.tex`, `docu/summary_of_changes_a70.tex` |
| 2013 (план) | — | «Python-версия процедур написания моделей (в работе, G. Tardini)» и «ASTRA как подпрограмма (долгосрочно)» — первое стало `pyparse/` в A8 | `docu/summary_of_changes_a70.tex`, раздел «Foreseen upgrades» |
| 2022-10-13 | 8.x | первый коммит в git MPCDF (`14208069`, G. Tardini, «Initial commit») | `git log` |
| 2023 | 8.x | FEQIS — решатель Града–Шафранова (E. Fable, 2023), второй после SPIDER; флаг flight simulator `-fs` (`bbfc1b25`, 2023-03-06) | `docu/FEQIS_manual.tex`, `git log` |
| февраль 2024 | 8.0 | «ASTRA-8.0 Addendum with changes and new features» (G. Tardini, E. Fable) | `docu/addendum_A8.tex` |
| июнь 2024 | 8.x | вывод в JSON с последующей склейкой в CDF вместо прежнего `A2CDF` (`02f7fae1` «Added json output!», `37617f1c`, `d5f0faee`) | `git log` |
| май 2026 | 8.x | заглушки для отсутствующих SPIDER, TGLF, QLKNN, RABBIT, TORBEAM, TORIC (`2ec99235`) | `git log` |
| июль–сентябрь 2026 | 8.x | связь с VMEC/STELLOPT, DKES/MONKES, переносимость на macOS и gfortran (`7698c2f7`, `b16ceeec`, `a4d11b69`, `16f08164`) | `git log` |
| январь 2026 | **8.6** | номер версии в баннере (`Version 8.6, January 2026`), выводится в заголовке окна | `repo_user_area/exe/version`, `src/graph/graph_utils.f90` (чтение `exe/version`) |

Статистика git (на `67281112`): 1806 коммитов с 2022-10-13; из них 1783 — Giovanni Tardini, единичные — Daniel Fajardo, Chuanren Wu, Michael Bergmann, Emiliano Fable и др. По годам: 2022 — 13, 2023 — 530, 2024 — 357, 2025 — 398, 2026 — 508. Тегов нет. История до 2022 года (ASTRA 7 и раньше) в этот репозиторий не входит.

## 3. Что изменилось в ASTRA 8 по сравнению с ASTRA 7

Источник — `docu/addendum_A8.tex` (введение и разделы 1–7) плюс сверка с кодом. Введение addendum формулирует главное так:

> Astra 8 is essentially a "new" code with respect to old Astra versions. The Code is now recast in F90, using modules and avoiding at all commons.

### 3.1. Сводная таблица

| область | ASTRA 7 | ASTRA 8 | где смотреть |
|---|---|---|---|
| язык ядра | Fortran 77 с `COMMON`-блоками (по addendum) | Fortran 90 модули, без `COMMON` | `src/astra/*.f90` |
| транслятор модели | собственный построитель (в этом репозитории его нет — не проверено, как он был устроен) | Python-парсер `pyparse/` → `src/tmp/*.f90` | `pyparse/parser_main.py`, [03-model-language.md](03-model-language.md) |
| структура установки | — | две области: репозиторная и пользовательская (`repo_user_area/` → копия в корень) | раздел 5 |
| сборка | — | `make` + `platform/env.<площадка>` + `exe/astra_rc` | [02-build-and-run.md](02-build-and-run.md) |
| запуск | `astra <эксперимент> <модель>` (по [`../../AGENTS.md`](../../AGENTS.md)) | `exe/as_exe -m <модель> -v <эксперимент> -s t0 -e t1 [-batch]`; без аргументов — повтор последнего расчёта | `repo_user_area/exe/as_exe.py` |
| синтаксис модели | A7 | 4 изменения (раздел 3.2) | `docu/addendum_A8.tex`, раздел 2 |
| равновесие | 3-моментное, SPIDER (в т.ч. свободная граница) | `IPEQL`: −2, −1, −6, 0, 1 (EMEQ), 4 (SPIDER, только заданная граница), **5 (FEQIS)**; свободная граница SPIDER больше недоступна | `docu/addendum_A8.tex`, раздел 3; [06-equilibrium.md](06-equilibrium.md) |
| описание машины для свободной границы | — | `exp/cnf/<dev>_description_in.json` → `exp/cnf/<dev>_description_full.json` (матрицы Грина) | `repo_user_area/exe/greenMatrices.py` |
| namelist-файлы модулей | — | `exp/nml/<эксперимент>` (или `exp/nml/<машина>`): RABBIT, NBI, TORBEAM, SPIDER, `equil_settings`, `fenix` | `docu/addendum_A8.tex`, `src/astra/read_input.f90::read_nml` |
| связь с внешней системой управления | — | режим flight simulator (`-fs`), интерфейс C. Wu (KIT); подпрограммы связи в репозиторий не входят | `docu/addendum_A8.tex`, раздел 5 |
| стелларатор | 1D-возможности с версии 3.0 | обобщённое уравнение диффузии тока, связь с VMEC/BOOZER/DKES (F. Solfronk, E. Fable) | `docu/astell_manual.tex`, `src/astra/a2vmec.f90`, `vmec/` |
| вывод | (не проверено для A7) | `ncdf_out/<exp><model>-<n>.json` по шагам → `ncdf_out/<exp><model>.CDF` | `src/astra/json_rw.f90`, `postProcess/json2cdf.py` |
| метаданные переменных | — | единый реестр `astra_variables.json` (читают и Fortran, и Python) | раздел 6.1 |

### 3.2. Изменения синтаксиса модели

Дословно по `docu/addendum_A8.tex`, раздел «Changes in the model file syntax» (подробный разбор — [03-model-language.md](03-model-language.md)):

1. Функции из `fnc/` задаются либо как массивы (`CARxx = QETOT`), либо как граничные интегралы (`CV1 = QETOTB`); вызывать их как функцию от `RHO` нельзя.
2. Значение массива в точке (например, `TE`, `TI`, `NE`) задаётся только через `AFX`: `TE(AFX(0.))` — на оси, `TE(AFX(0.5))` — на половине малого радиуса.
3. Если задано начальное условие (`NE = NEX`), граничное условие ставится автоматически как `NEB = NEX(граница из команды EQ[])`; явное `NEB` имеет приоритет.
4. Формулу из `fml/` и команду `GRAD` нельзя использовать в одной строке.

Отдельно: в репозитории лежат модели в формате A7 (`repo_user_area/equ/showdata`, `repo_user_area/exp/readme`), которые в A8 не запускаются — см. [`../../installer/docs/DETAILED.md`](../../installer/docs/DETAILED.md), раздел 7, и раздел 10 ниже.

### 3.3. Режимы равновесия (`IPEQL`)

`IPEQL` задаётся в блоке параметров файла `equ/log/<модель>`. По `docu/addendum_A8.tex`, раздел «Equilibrium solvers»:

| `IPEQL` | смысл |
|---|---|
| −2 | цилиндрическое аналитическое равновесие |
| −1 | метрические коэффициенты из файла эксперимента (`SHX`, `SHIVX`, `SQUAX`, `ELX`, `TRX`, `G33X`, `IPOLX`, `VRX`, `SLATX`, `G11X`, `G22X`, `DRODAX`) |
| −6 | то же, но из файла `input_metric.dat` в корне установки, уже на сетке ASTRA |
| 0 | фиксированное начальное приближение, без внешнего равновесия |
| 1 | 3-моментное равновесие EMEQ (`src/astra/emeq.f90`) |
| 2, 3 | не используются |
| 4 | SPIDER, заданная граница (закрытый код IPP) |
| 5 | FEQIS, заданная или свободная граница (`src/feqis/`) |

При `flightsim` режимы −1 и −6 обходятся: величины должны прийти извне. На коммите `e66f6ba1` (2026-10-01) все тестовые модели, кроме `equ/flux_spider`, переведены на `IPEQL=5` (сообщение коммита «Moved test cases to IPEQL=5»). Подробно — [06-equilibrium.md](06-equilibrium.md).

## 4. Авторы, сопровождающие, лицензия

**Создатели**: Г. В. Переверзев (G. V. Pereverzev) и П. Н. Юшманов (P. N. Yushmanov) — авторы ASTRA и руководства IPP 5/98. В баннере `repo_user_area/exe/version` (выводится программой) — «Creators: G. V. Pereverzev and Y. P. Yushmanov», «in memoriam of Prof. G. V. Pereverzev».

**Разработчики ASTRA 8** (`repo_user_area/exe/version`, `README.md`):

- Giovanni Tardini (`giovanni.tardini@ipp.mpg.de`) — сопровождающий, автор подавляющего большинства коммитов, Python-парсера, сборочной системы;
- Emiliano Fable (`emiliano.fable@ipp.mpg.de`) — сопровождающий, автор ASTRA 7.0 addendum, FEQIS, стеллараторной части (вместе с E. Buglione-Ceresa — `vmec/python/run_dkes_parallel.py`);
- дополнительные благодарности: A. Merle (CRPP), A. Polevoi (ITER Organization);
- стеллараторное уравнение диффузии тока — Fabian Solfronk (`docu/astell_manual.tex`);
- интерфейс flight simulator — C. Wu, KIT (`docu/addendum_A8.tex`).

**Организация**: Max-Planck-Institut für Plasmaphysik (IPP), Garching. Оригинальный репозиторий — <https://gitlab.mpcdf.mpg.de/git/astra.git> (анонимный HTTPS-доступ к нему есть; закрытые модули IPP — только по SSH для аккаунтов IPP).

**Лицензия**: GNU LGPL v2.1 или новее — файл `LICENSE`, `README.md`. Заголовок исходников, например `src/astra/astra_main.F90`, строки 3–8: «Copyright (C) 2026 Institut fuer Plasmaphysik ... GNU Lesser General Public License ... version 2.1». Правила обращения с форком (сохранять `LICENSE` и копирайты, помечать модифицированную версию) — в [`../../docs/local-fork.md`](../../docs/local-fork.md).

Внешние модули имеют **свои** лицензии и в основном не входят в репозиторий: RABBIT, TORBEAM, SPIDER, STRAHL — закрытые коды IPP; TGLF/NEO (GACODE), QuaLiKiZ, QLKNN, json-fortran — открытые. См. [07-heating-and-external-modules.md](07-heating-and-external-modules.md).

## 5. Репозиторная и пользовательская области

### 5.1. Идея

`docu/addendum_A8.tex`, раздел 1: после клонирования код имеет **два слоя**:

1. корень (`a8code/`) — исходники ядра (`src/`), парсер, скрипты;
2. `a8code/repo_user_area/` — рабочие каталоги с готовыми формулами, подпрограммами, примерами моделей и данных.

`install.sh` копирует **все** подкаталоги `repo_user_area/` в корень. Пользователь работает уже с копиями (`equ/`, `exp/`, `sbr/`, ...), может добавлять свои файлы и менять поставляемые, не трогая репозиторный слой. Копии игнорируются git (`.gitignore`), поэтому `git pull` не конфликтует с работой пользователя. Обратная сторона: обновления примеров и сборочных скриптов из `repo_user_area/` **не попадают** в рабочие копии сами — их нужно копировать заново (`install.sh` или установщик), а `install.sh` без `-safe` перезаписывает одноимённые файлы (подробно — [02-build-and-run.md](02-build-and-run.md), раздел 2).

### 5.2. Что копируется

| каталог в `repo_user_area/` | записей | назначение после копирования в корень | подробно |
|---|---|---|---|
| `equ/` | 27 моделей + `log/` (26 файлов) | модели (язык ASTRA) и их параметры `equ/log/<модель>` | [03](03-model-language.md), [09](09-user-area-catalog.md) |
| `exp/` | 6 файлов + `cnf/`, `icrh/`, `nbi/`, `nml/` | файлы эксперимента, описания машин, namelist'ы модулей | [04](04-experiment-data.md) |
| `udb/` | 3 разряда (`33040_2500`, `34954`, `36982_3400`) | U-файлы и файлы границы для примеров | [04](04-experiment-data.md) |
| `fml/` | 311 | формулы (фрагменты Fortran, подставляемые парсером) | [09](09-user-area-catalog.md) |
| `fnc/` | 117 | функции `*.f90` → `lib/libfnc.a` | [09](09-user-area-catalog.md) |
| `sbr/` | 37 | подпрограммы `*.f90/*.F90` → `lib/libsbr.a` (в т.ч. интерфейсы RABBIT, TORBEAM, TGLF, QLK, NEO, STRAHL) | [07](07-heating-and-external-modules.md), [09](09-user-area-catalog.md) |
| `xpr/` | 4 | исходники IPC-программ (`tglf_interf.f90`, `qlk_interf.f90`, `neo_interf.f90`, `tra_interf.c`) | [07](07-heating-and-external-modules.md) |
| `exe/` | 14 | скрипты сборки и запуска | [02](02-build-and-run.md) |
| `preProcess/` | 2 | `sf2astra.py`, `imas2astra.py` — подготовка данных из shotfile AUG / IMAS (IPP-специфично) | [04](04-experiment-data.md) |
| `pyparse/` | 1 | `iondens_ass.py` — пользовательский фрагмент присваивания `NI`, кладётся **внутрь** репозиторного `pyparse/` | раздел 6.3 |

Примечание про `equ/`: у `showdata` нет файла `equ/log/showdata`, поэтому эта модель упадёт ещё на чтении параметров (`src/astra/read_input.f90::readInput`, сообщение `>>> Error: file "equ/log/showdata" missing`), независимо от проблем формата A7.

### 5.3. `.gitignore`

Корневой `.gitignore` (весь файл, 2026-10-01) игнорирует:

- пользовательскую область: `/equ`, `/exp`, `/udb`, `/exe`, `/fml`, `/fnc`, `/sbr`, `/xpr`, `/preProcess`, `pyparse/iondens_ass.py`;
- результаты сборки: `mod`, `obj`, `/lib`, `/bin`, `/src/tmp`, `**/*genmod.f90`, `*.pyc`, `__pycache__`;
- результаты расчёта и рабочие каталоги: `/ncdf_out`, `/dat`, `/tmp`, `/ipc`, `/vmec/dat`, `/tglf`, `/qualikiz`, `/neo`, `/toric`, `/helena`, `/sicas`, `profiles.dat`, `lastruntime.qlk`, `out.neo.*`, `regressions/*.log`, `preProcess/memo`, `slurm-*`;
- прочее: `user`, `*_backup` (резервные копии `install.sh -safe`), `**/*~`.

**Находка**: в `.gitignore` остались маркеры неразрешённого слияния:

```
<<<<<<< HEAD
/sicas
=======
>>>>>>> stars3d
```

git воспринимает эти строки как обычные шаблоны (файлы с такими именами вряд ли появятся), так что на работу это не влияет, но говорит о небрежном слиянии ветки `stars3d`.

## 6. Карта дерева исходников

Число отслеживаемых git файлов по каталогам верхнего уровня (`git ls-files`, `67281112`): `repo_user_area` — 592, `strahl` — 153, `src` — 56, `docu` — 29, `rabbit` — 29, `platform` — 22, `pyparse` — 13, `regressions` — 12, `vmec` — 10, `postProcess` — 8, `inc` — 1, `python` — 1, плюс 14 файлов в корне. Размеры на диске: `regressions` 104 МБ, `rabbit` 13 МБ, `strahl` 8.2 МБ, `repo_user_area` 3.2 МБ, `docu` 1.8 МБ, `src` 1.3 МБ.

### 6.1. Файлы в корне

| файл | назначение | подробно |
|---|---|---|
| `README.md` | клонирование, `./install.sh`, `exe/as_exe`, список поддерживаемых площадок, сопровождающие, лицензия | — |
| `LICENSE` | текст GNU LGPL v2.1 | — |
| `.gitignore` | см. раздел 5.3 | — |
| `install.sh` | первичная установка: копирование `repo_user_area/*`, `chmod`, `make clean`, тестовый расчёт `flux_feqis`/`aug34954`; режим `-safe` с вопросами | [02](02-build-and-run.md), разд. 2 |
| `get_platform` | печатает имя площадки по пути установки, `uname -n` и ОС; иначе `mypc`; суффикс `_${ASTRA_COMPILER}` | [02](02-build-and-run.md), разд. 3; [`platforms.md`](../../docs/astra/platforms.md) |
| `clean.sh` | `make -f exe/Makefile clean` из корня установки | [02](02-build-and-run.md), разд. 9 |
| `compReg.sh` | прогон набора регрессионных расчётов и сравнение с `regressions/*.CDF` | [08](08-graphics-and-output.md) |
| `compareRegressions.py` | сравнение `ncdf_out/<exp><model>.CDF` с `regressions/<exp><model>.CDF`; допуск по умолчанию `1e-7` (функция `compare`); лог в `regressions/<exp><model>.log` | [08](08-graphics-and-output.md) |
| `plot_compare.py` | набросок: рисует `pprime` на последнем шаге для `aug34954fluxes` (эталон против нового); имя файла зашито, модели `fluxes` в `repo_user_area/equ/` нет | — |
| `modules_install.sh` | интерактивная сборка внешних модулей (Intel oneAPI, CMake, json-fortran, NetCDF, RABBIT, TORBEAM, SPIDER, QuaLiKiZ, QLKNN, TGLF/NEO, STRAHL) в `$ASTRA_EXT` для кластеров | [02](02-build-and-run.md), разд. 6 |
| `modules_inst_loc.sh` | то же для «локальных» машин (ветки Linux/apt и macOS/Homebrew), плоская раскладка `$ASTRA_EXT` | [02](02-build-and-run.md), разд. 6 |
| `requirements.txt` | человекочитаемый список требований (не формат pip) | [02](02-build-and-run.md), разд. 7 |
| `astra_variables.json` | реестр всех переменных ASTRA с единицами и описаниями (см. ниже) | [05](05-solver-core.md) |
| `memo_stella` | памятка: команда запуска стеллараторного примера и нужные подкаталоги STELLOPT | ниже |

**`astra_variables.json`** — словарь из 14 групп (число записей): `variables` (135), `variables_x` (135), `constants` (162), `control` (15), `internInt` (33), `internDbl` (46), `profiles` (454), `profiles_x` (106), `strahl` (13), `equil_signals` (8), `equil_profiles` (28), `equil_rect` (4), `equil_rz2d` (5), `equil_coord` (6). Каждая запись — `{"units": ..., "desc": ...}`. Читают его:

- Fortran при старте — `src/astra/astra_main.F90:45` (`call read_metadata()`, реализация в `src/astra/json_vars.f90`, имя файла в строке 53) и `associate_pointers()` из сгенерированного `src/tmp/associate_pointers.f90`;
- парсер — `pyparse/equ_parser.py:48`, `pyparse/exp_parser.py:99`, `pyparse/code_gen.py`;
- конвертер вывода — `postProcess/json2cdf.py:35` (метаданные для CDF).

Значит, добавить новую «встроенную» переменную — это в первую очередь правка `astra_variables.json` (подробно — [05-solver-core.md](05-solver-core.md)). Файл читается по относительному пути, поэтому исполняемый файл должен запускаться из корня установки (это делает `exe/Build`).

**`memo_stella`** (весь файл):

```
exe/as_exe -m transpstell -v steltest -s 0. -e 4. -tpause 0. -nodes 8

Required STELLOPT subdirs:
    BOOZ_XFORM
    VMEC2000
    DKES
```

### 6.2. `src/` — ядро (Fortran/C)

| каталог | файлов | содержимое | библиотека | подробно |
|---|---|---|---|---|
| `src/astra/` | 28 | главная программа `astra_main.F90`, шаг по времени `stepup.F90`, транспортный решатель `transport_solver.f90`, чтение входа `read_input.f90`, метрика `metrics.f90`, EMEQ `emeq.f90`, связь с равновесием `gs_solver.F90`, VMEC/DKES `a2vmec.f90`, `a2dkes.f90`, `vmec_indata.f90`, `stella_module.f90`, JSON-ввод/вывод `json_rw.f90`, `json_vars.f90`, IPC `ipc_control.c`, описание машины `machine_config.f90`, скаляры/статус/функции | `lib/libastra.a` (кроме `astra_main.F90` и `stepup.F90`, которые компилируются при линковке) | [05](05-solver-core.md), [06](06-equilibrium.md) |
| `src/feqis/` | 13 | FEQIS: заданная (PBE) и свободная (FBE) граница, цепи катушек `circuit.f90`, связь с транспортом `transport2fbe.f90` | `lib/libfeqis.a` | [06](06-equilibrium.md) |
| `src/graph/` | 4 | X11-окно: `Astra2XW.c`, `dialogx.c`, `graph_utils.f90`, `gui_interaction.f90` | `lib/libgraph.a` (линкуется только в GUI-сборке) | [08](08-graphics-and-output.md) |
| `src/misc/` | 5 | общие утилиты: строки, числа, `pi_const`, `numerical_tools`, `imas_ids` (структуры IMAS) | `lib/libmisc.a` | [05](05-solver-core.md) |
| `src/nbi/` | 6 | встроенный NBI-модуль (сечения, инжекция) + `ReadMe.txt` | `lib/libnbi.a` | [07](07-heating-and-external-modules.md) |

Кроме того, при сборке появляется `src/tmp/` (игнорируется git) — сгенерированный парсером Fortran-код текущей модели (см. раздел 8 и [02](02-build-and-run.md)).

### 6.3. `pyparse/` — транслятор моделей

13 файлов Python. Точка входа — `pyparse/parser_main.py` (`-equ <путь к модели> -exp <путь к эксперименту>`), вызывается из `exe/Makefile` (переменная `PY_TMP_GEN`, строка 173) и пишет в `./src/tmp`:

| модуль | роль |
|---|---|
| `parser_main.py` | связывает разбор модели и эксперимента, вычисляет `nr_x_max` (максимальное число радиальных точек в профилях эксперимента), вызывает генератор и пишет файлы |
| `equ_parser.py` | разбор файла модели `equ/<m>` (по реестру `astra_variables.json`, `fml/`, `fnc/`) |
| `exp_parser.py` | разбор файла эксперимента `exp/<e>`, ссылки на U-файлы |
| `ufiles.py` | чтение U-файлов |
| `parse_as.py` | лексика языка ASTRA, форматирование Fortran (`write_fortran`) |
| `eqns.py`, `eqns_init.py` | генерация кода транспортных уравнений и инициализации |
| `code_gen.py` | сборка итоговых подпрограмм (`CODE_GEN`) |
| `const_text.py` | шаблоны неизменяемых частей Fortran-кода |
| `fml.py` | подстановка формул из `fml/`; есть отладочный CLI (`-s`, `-f`) |
| `config.py` | константы (`NCVA=512`, `NFML=500`, `NSBMX=20`, `NRW=128`), список уравнений `eqn_list`, имена коэффициентов, потоков и граничных условий для каждого уравнения |
| `rw_for.py` | чтение/запись Fortran-форматированных чисел (дублирует `postProcess/rw_for.py`) |
| `test_parse.py` | отладочный разбор одного оператора |
| `iondens_ass.py` | **не в git** (`.gitignore`), копируется из `repo_user_area/pyparse/`; задаёт присваивание `NI = F1 + ... + F9` |

Подробно — [03-model-language.md](03-model-language.md).

### 6.4. `postProcess/` и `python/`

| файл | назначение |
|---|---|
| `postProcess/json2cdf.py` | склейка `ncdf_out/<exp><model>-<n>.json` в `ncdf_out/<exp><model>.CDF` (NetCDF-3 через `scipy.io.netcdf_file`); вызывается из `exe/as_exe.py` (`json_concat`) и отдельно (`-e <exp><equ>`) |
| `postProcess/json2eqdsk.py` | JSON-вывод → gEQDSK (`-e`, `-t` — индекс времени) |
| `postProcess/eqdsk.py` | чтение/запись gEQDSK, таблица COCOS |
| `postProcess/acdf2imas.py` | CDF ASTRA → IMAS IDS (нужен модуль `imas`; бэкенды HDF5, MDSPLUS, ASCII) |
| `postProcess/plotdkesrun.py` | графики результатов DKES |
| `postProcess/readxmxn.py` | чтение мод `xm`, `xn` из файла VMEC |
| `postProcess/rw_for.py` | Fortran-формат чисел |
| `postProcess/test_pressure.py` | отладочный скрипт (назначение не изучалось) |
| `python/eqdsk_to_TORIC_gs.py` | gEQDSK → файл равновесия TORIC (контуры, Фурье-разложение, сглаживание, уточнение Града–Шафранова) |

Подробно — [08-graphics-and-output.md](08-graphics-and-output.md).

### 6.5. Данные внешних модулей: `rabbit/`, `strahl/`, `vmec/`

| каталог | содержимое | подробно |
|---|---|---|
| `rabbit/limiters/` | контуры лимитеров для RABBIT: `limiter_{aug,d3d,dtt,jet,tcv,vns,COMU,future}.dat` | [07](07-heating-and-external-modules.md) |
| `rabbit/tables_ITERDEMO/`, `rabbit/tables_highRes/` | таблицы сечений (`dd/*_table.dat`, `einstein.dat`) — данные для RABBIT; сам код RABBIT закрыт и не входит | [07](07-heating-and-external-modules.md) |
| `strahl/` | атомные данные STRAHL `*.atomdat` (Ar, B, Be, C, Cu, Fe, H, He, Kr, Li, N, Ne, Ni, O, W, Xe), `atomdat/`, `strahl.control`, `pec_files` (пути к ADAS вида `/u/rld/adas` — IPP-специфично), `pstrahl` (C-shell-скрипт запуска IDL) | [07](07-heating-and-external-modules.md) |
| `vmec/python/` | `vmec2bin.py` (самодостаточная обработка VMEC без STELLOPT), `vmec.py`, `extract_boozer_data.py`, `parse_fortran_nml.py`, `pellet.py`, `run_dkes_parallel.py`, `run_monkes_parallel.py` | [06](06-equilibrium.md) |
| `vmec/templates/` | `vmecinput_template.dat`, `stell_files.nml`, `true_surf.txt` | [06](06-equilibrium.md) |

### 6.6. `platform/`, `inc/`, `regressions/`, `docu/`

| каталог | содержимое | подробно |
|---|---|---|
| `platform/` | 22 файла `env.<площадка>`: `aug`, `columbia`, `cz`, `darwin`, `docker`, `freia`, `gway`, `hgw`, `iter`, `lac`, `mit`, `ncku`, `omega`, `perlmutter`, `puhti`, `rat2`, `sevilla`, `tohtori`, `tok`, `tok_gcc`, `tok_ifx`, `w7x`. Каждый задаёт компилятор, флаги, MKL, `ASTRA_EXT`, иногда SLURM-очередь. Обычного Linux-ПК (`env.mypc`) нет | [`platforms.md`](../../docs/astra/platforms.md) |
| `inc/Astra.h` | 22 строки: заголовки System V IPC (семафоры, разделяемая память), `union semun` для не-macOS; подключается при компиляции C-файлов (`src/astra/ipc_control.c`, `xpr/tra_interf.c`) | [07](07-heating-and-external-modules.md) |
| `regressions/` | 12 эталонных CDF (Intel+MKL, кластер IPP): `AUG33040_2500fbe`, `AUG36982_3400tglf_pid`, `aug34954_tflux_feqis`, `aug34954astrahl_simple_msp`, `aug34954flux_{cuas,feqis,nbi,neo_tglf,spider}`, `aug34954{qlk,qlknn,tglf}` | [08](08-graphics-and-output.md); [`DETAILED.md`](../../installer/docs/DETAILED.md), разд. 5.3 |
| `docu/` | LaTeX-исходники и PDF документации | раздел 7 |

## 7. Документация в `docu/`

`docu/compilingmain` — рецепт сборки: `latex` дважды, `dvips`, `ps2pdf` для четырёх документов: `main.tex`, `addendum_EF.tex`, `addendum_A8.tex`, `runeq.tex`. Готовые PDF лежат рядом.

| файл | строк | что это | входит в |
|---|---|---|---|
| `main.tex` → `main.pdf` | 96 | основное руководство: G. V. Pereverzev, P. N. Yushmanov, «ASTRA — Automated System for TRansport Analysis», IPP 5/98, февраль 2002 (версия 5.3) | — |
| `section12.tex` | 722 | разделы 1–2: введение, история, обзор кода, диаграммы «пользователь ↔ ASTRA» и «поток сборки» | `main.tex` |
| `section3.tex` | 1561 | раздел 3: физическая база — координаты и усреднение, ток и поле, равновесие, электрическое поле, омический нагрев, уравнения переноса, источники, начальные и граничные условия, замыкание, вспомогательные уравнения, нейтралы | `main.tex` |
| `section4.tex` | 3208 | раздел 4 (Reference guide): единицы и переменные, уравнения в обозначениях ASTRA, связь транспорта и равновесия, переменные, выражения, встроенные функции, подключаемые подпрограммы, интерфейсы нагрева | `main.tex` |
| `section5.tex` | 3440 | раздел 5 (User's guide): концепция, задание модели, задание данных, формат файлов данных, краткое руководство по работе | `main.tex` |
| `bib.tex` | 57 | библиография основного руководства | `main.tex` |
| `definitions.tex` | 79 | макросы LaTeX | `main.tex`, `addendum_*.tex` |
| `addendum_EF.tex` → `addendum_EF.pdf` | 570 | E. Fable, «ASTRA-7.0 Addendum with changes and new features» (2012/2013) | — |
| `norm_grid_add_EF.tex` | 181 | нормированная сетка как опорная | `addendum_EF.tex` |
| `mom_transp_eq_add_EF.tex` | 312 | уравнение тороидального момента | `addendum_EF.tex` |
| `coupling_with_2d_equil_add_EF.tex` | 362 | связь с 2D-равновесием | `addendum_EF.tex` |
| `icp_add_EF.tex` | 76 | параллельные процессы IPC | `addendum_EF.tex` |
| `stiff_transport_EF.tex` | 76 | численная схема для жёсткого транспорта | `addendum_EF.tex` |
| `newton_solver.tex` | 152 | Ньютоновская схема для стационара | `addendum_EF.tex` |
| `additional_packages_add_EF.tex` | 135 | дополнительные пакеты ASTRA 7.0 | `addendum_EF.tex` |
| `additional_changes_add_EF.tex` | 58 | прочие изменения 7.0 | `addendum_EF.tex` |
| `summary_of_changes_a70.tex` | 59 | сводка изменений 7.0 относительно 6.2.1 и планы | `addendum_EF.tex` |
| `addendum_A8.tex` → `addendum_A8.pdf` | 697 | G. Tardini, E. Fable, «ASTRA-8.0 Addendum with changes and new features», февраль 2024: двухслойная структура, CLI, синтаксис, `IPEQL`, генерация описания машины для FEQIS, namelist'ы, связь с FEQIS (полный список параметров `transport2fbe`), flight simulator | — |
| `FEQIS_manual.tex` | 340 | FEQIS: входные файлы (`exp/equ/<машина>/coil.dat`, ...) — **но дальше, со строки ~30, идёт текст про уравнение момента из A7** (см. раздел 10) | `addendum_A8.tex` |
| `astell_manual.tex` | 234 | стелларатор: обобщённое уравнение диффузии тока, \(j_\parallel\) в неаксиальной геометрии, перенос частиц и энергии, новые переменные, равновесие, \(\rho_B\) (F. Solfronk) | `addendum_A8.tex` |
| `Stellarator_Astra.tex` | 340 | несмотря на имя — почти точная копия `FEQIS_manual.tex` (126 строк `diff`, в основном пробелы в концах строк); про стелларатор не написано | ни во что не включён |
| `IPC.tex` | 563 | отдельный документ о параллельных процессах ASTRA 6.1 (`NEUT & : 0.1 : 2.2 : 3.3;`); подключает `../input/colors`, `../input/math`, которых в репозитории нет, — в таком виде не собирается | ни во что не включён |
| `runeq.tex` → `runeq.pdf` | 80 | «Generic transport equation»: общий вид уравнения переноса ASTRA с конвективными членами сетки, явная и неявная формы | — |
| `def.tex` | 64 | макросы для `runeq.tex` | `runeq.tex` |

Документация описывает идеи и физику хорошо, но в деталях реализации A8 часто отстаёт от кода (раздел 10). Описание переменных, уравнений и численной схемы — в [05-solver-core.md](05-solver-core.md).

## 8. Поток данных: от модели до CDF

### 8.1. Схема

```
          ┌───────────────── пользователь ─────────────────┐
          │                                                 │
   equ/<m>           equ/log/<m>        exp/<e>     exp/nml/<e|dev>   exp/cnf/<dev>_description_in.json
   (модель)          (параметры)        (данные)    (namelist'ы)      (машина, катушки)
      │                   │                │              │                    │
      │                   │                ├─ udb/<shot>/* (U-файлы)           │
      │                   │                │                                   ▼
      ▼                   │                ▼                       exe/greenMatrices.py
  ┌────────────────────────────────────────────┐                  (при каждом as_exe, если
  │ pyparse/parser_main.py -equ ... -exp ...   │                   вход новее выхода)
  │  + astra_variables.json, fml/*, fnc/*      │                           │
  └────────────────────────────────────────────┘                           ▼
      │                                                exp/cnf/<dev>_description_full.json
      ▼                                                (матрицы Грина для FEQIS)
  src/tmp/*.f90  (associate_pointers, inivar, detvar, postep,        │
                  ininam, init_converge_step, eqns_inc,              │
                  [astra_out, set_graph_names — только GUI])         │
      │                                                              │
      ▼                                                              │
  gfortran/ifx: src/astra/astra_main.F90 + stepup.F90 (или -workflow)│
               + src/tmp/*.f90                                       │
               + lib{sbr,nbi,graph,astra,fnc,feqis,misc}.a           │
               + json-fortran, MKL/OpenBLAS, [X11], [внешние модули] │
      │                                                              │
      ▼                                                              │
  bin/<m>_gui.exe  или  bin/<m>_batch.exe                            │
      │   запуск из корня установки (exe/Build делает cd $AWD)       │
      │   читает: tmp/astra.nml, astra_variables.json, equ/log/<m>,  │
      │           exp/nml/..., exp/<e>, udb/..., exp/cnf/... ◄───────┘
      ▼
  ncdf_out/<e><m>-1.json, -2.json, ...   (каждые DPOUT секунд модельного времени)
      │
      ▼  postProcess/json2cdf.py (вызывается из exe/as_exe.py)
  ncdf_out/<e><m>.CDF
      │
      ├─► compareRegressions.py (сравнение с regressions/*.CDF)
      ├─► postProcess/json2eqdsk.py, acdf2imas.py
      └─► внешние просмотрщики CDF
```

### 8.2. Пояснения

- **Модель** (`equ/<m>`) определяет, *какой* код будет собран. Исполняемый файл называется по модели: `bin/<m>_gui.exe` / `bin/<m>_batch.exe` (`repo_user_area/exe/Build`, строки 35–42). Эксперимент в имя не входит.
- **Эксперимент** (`exp/<e>`) влияет на сборку только числом `nr_x_max` (размер массивов экспериментальных профилей; `pyparse/parser_main.py:30`, `pyparse/code_gen.py:98`); остальное читается во время выполнения (`src/astra/read_input.f90::read_exp`). Следствие для make — в [02](02-build-and-run.md), раздел 10.
- **Параметры** (`equ/log/<m>`) — значения переменных (`AB`, `BTOR`, `IPL`, ...), констант и управляющих величин; читаются во время выполнения (`src/astra/read_input.f90`, строки 84–108). Их можно менять без пересборки.
- **Командная строка** записывается в `tmp/astra.nml` (namelist `&astra_log`) и читается программой (`src/astra/read_input.f90::read_nml`); время начала/конца из командной строки перекрывает `equ/log`.
- **Путь установки** вшивается в исполняемый файл: сгенерированный `src/tmp/ininam.f90` содержит строку `AWD = "<путь>"` (видно в `~/astra/a8/src/tmp/ininam.f90:12`), по которому пишутся JSON (`src/astra/json_rw.f90:413`). Перенос каталога установки требует пересборки.
- **Вывод**: каждые `DPOUT` секунд модельного времени `write_ajson` пишет очередной `ncdf_out/<e><m>-<n>.json` (`src/astra/astra_main.F90`, строки 122–129); после выхода программы `exe/as_exe.py` склеивает их в CDF (`postProcess/json2cdf.py::json_concat`). Пример: `flux_feqis`/`aug34954`, t = 4–5 с, 21 шаг по 50 мс — 6 JSON-файлов, CDF с `nt=6, nrho=91` (`~/astra/test_flux_feqis.log`).

## 9. Что вычисляет один расчёт

Здесь — только концептуальная картина; уравнения, переменные и численная схема — в [05-solver-core.md](05-solver-core.md), равновесие — в [06-equilibrium.md](06-equilibrium.md).

### 9.1. Набор уравнений

ASTRA решает систему одномерных (по радиальной координате \(\rho\), в A7+ — по нормированной \(x=\sqrt{\Phi/\Phi_b}\)) диффузионно-конвективных уравнений, усреднённых по магнитным поверхностям. Список уравнений, которые может включить модель, — `pyparse/config.py`, `eqn_list`:

| уравнение | величина | коэффициенты (из `config.py::coeff_d`) | источники (`flux_d`) | граница (`bnd_d`) |
|---|---|---|---|---|
| `NE` | плотность электронов | `DN`, `HN`, `XN`, `CN`, `DVN`, `DSN` | `SN`, `SNN` | `NEB`, `QNB`, `QNNB` |
| `TE` | температура электронов | `DE`, `HE`, `XE`, `CE`, `DVE`, `DSE` | `PE`, `PET` | `TEB`, `QEB`, `QETB` |
| `TI` | температура ионов | `DI`, `HI`, `XI`, `CI`, `DVI`, `DSI` | `PI`, `PIT` | `TIB`, `QIB`, `QITB` |
| `CU` | ток / полоидальный поток (\(\mu = 1/q\)) | `DC`, `HC`, `XC`, `CC`, `CD` | `MV`, `MU` | `IPL`, `MUB` |
| `UPAR` | параллельная скорость (тороидальный момент) | `CNPAR`, `XUPAR`, `RUPAR`, `RUPYR`, `RUPFR` | `TTRQ` | `UPARB`, `TTRQB` |
| `F0` … `F9` | 10 вспомогательных величин (примеси, нейтралы и т.п.) | `DFj`, `VFj`, `GFj`, `DVFj`, `DSFj` | `SFj`, `SFFj` | `FjB`, `QFjB`, `QFFjB` |
| `Equil` | равновесие (выбор решателя через `IPEQL`) | — | — | — |

Каждое уравнение в модели может быть «решаемым» (эволюционирует по уравнению переноса) или «заданным» (профиль присваивается, например, из эксперимента: `TE = TEX`). Общий вид уравнения — `docu/runeq.tex`, физика — `docu/section3.tex`, обозначения ASTRA — `docu/section4.tex`.

### 9.2. Ход расчёта

По `src/astra/astra_main.F90`:

1. **Инициализация** (строки 42–53): чтение реестра переменных (`read_metadata`), связывание указателей (`associate_pointers`, сгенерировано), значения по умолчанию (`scalars_init`), имена из модели (`ininam`), чтение входа (`readInput`: `tmp/astra.nml`, `exp/nml`, `exp/nbi`, `equ/log/<m>`, `exp/<e>` с U-файлами), описание машины (`config_read`).
2. **Окно** (строки 55–61): только в GUI-сборке (`#ifdef X11`) и если `TASK /= "BGD"`.
3. **Начальное состояние** (строки 63–108): либо restart из JSON-файла (`-re N` → `read_ajson`), либо итерации до самосогласования: `INIVAR` (начальные профили из модели), `DETVAR`, `eqguess`, затем цикл `set_x_scalars` → `DETVAR` → `DEFARR` → `INIT_CONVERGE_STEP` → `METRIC` до сходимости `IFTREQ(ATREQ)`. При `IPEQL == 5` инициализируется FEQIS (`transport2fbe_init`).
4. **Цикл по времени** (строки 122–129): пока `TIME < TEND`, каждые `DPOUT` — запись JSON, затем `STEPUP()` — один шаг: коэффициенты переноса и источники из модели, решение уравнений, обновление равновесия, вызовы подпрограмм `sbr/` с их периодичностью, графика. `STEPUP` можно заменить своим файлом (`-workflow`, см. [02](02-build-and-run.md)).
5. **Выход** (строка 131): `CPU_report('>>> ASTRA normal exit >>>')` — отчёт о времени по подпрограммам. Именно эта строка, а не код возврата, — признак успешного расчёта.

### 9.3. Предсказательный и интерпретационный режимы

Это не флаг, а свойство модели: если `TE`, `TI`, `NE` решаются по уравнениям с модельными коэффициентами переноса — расчёт предсказательный; если профили заданы из эксперимента, а вычисляются, например, только ток и равновесие (как в тестовой `flux_feqis`, где `TE`, `TI`, `NE` берутся из U-файлов) — интерпретационный. Руководство (`docu/section12.tex`) подчёркивает, что ASTRA строит коды «for predictive or interpretative transport modeling, for stability analysis or for processing experimental data».

## 10. Расхождения документации и кода

Найдено при сверке `docu/addendum_A8.tex` и других файлов с кодом на `67281112`. Код первичен.

| документация | факт в коде | где |
|---|---|---|
| в `src/` есть `for`, `feqis`, `spider`, `nbi` и каталог `main/` с `astra_main.f90`, `stepup.f90` и списками переменных | `src/` = `astra`, `feqis`, `graph`, `misc`, `nbi`; главная программа — `src/astra/astra_main.F90`, шаг — `src/astra/stepup.F90`; списки переменных — `astra_variables.json` в корне | `docu/addendum_A8.tex`, разд. 1 |
| CLI: `-m -v -s -e -dev -tpause -f -b -i <template>` | флаги называются `-fs` (flight simulator) и `-batch`; `-f` и `-b` работают только как сокращения argparse (проверено на Python 3.12 с тем же набором опций: `-b` → `-batch`, `-f` → `-fs`; а вот `-d` неоднозначен: `-dev`/`-debug`); опции `-i` (INCLUDE-шаблоны) нет | `repo_user_area/exe/as_exe.py`, строки 46–68 |
| «при отсутствии библиотек закомментировать вызовы xpr в `exe/Build`» | IPC-программы собираются в `exe/Makexpr` (цель `all`), `Build` только вызывает его | `repo_user_area/exe/Build`, строки 51–58 |
| описание машины для FEQIS: `exp/cnf/machine_description_out.<dev>`, генерировать вручную `green/greenMatrices.py` | файл `exp/cnf/<dev>_description_full.json`; генератор `exe/greenMatrices.py` вызывается автоматически из `exe/as_exe.py` (`greenMatrices.main()`, строка 43) для всех `exp/cnf/*_description_in.json`; FEQIS читает его в `src/feqis/feqis_solvers.f90:561` | `docu/addendum_A8.tex`, разд. 4 |
| пример связи с FEQIS — `repo_user_area/sbr/eqctst.f90` | такого файла нет | `docu/addendum_A8.tex`, разд. «Coupling to FEQIS» |
| входные файлы FEQIS в `exp/equ/<машина>/` | в `repo_user_area/exp/` нет каталога `equ/`; описание машины берётся из `exp/cnf/` (JSON); старый формат, по-видимому, устарел (не проверено) | `docu/FEQIS_manual.tex` |
| `FEQIS_manual.tex` про FEQIS | после ~30 строк про FEQIS идёт копия раздела «Toroidal momentum transport equation» ASTRA 7 (`\subsection{Derivation of the equation implemented in Astra--7.0}` и т.д.) — в `addendum_A8.pdf` это попадает в главу «The new equilibrium solver FEQIS» | `docu/FEQIS_manual.tex`, строки 32–340 |
| `Stellarator_Astra.tex` | дубликат `FEQIS_manual.tex`, нигде не подключён; стелларатор описан в `astell_manual.tex` | `docu/` |
| `astell_manual.tex`: уравнения решаются в `src/for/runeq.f90` | каталога `src/for/` нет; решатель — `src/astra/transport_solver.f90` (предположительно; см. [05](05-solver-core.md)) | `docu/astell_manual.tex:66` |
| addendum: «аргументы необязательны, используются последние значения» | верно только если уже есть `tmp/astra.nml`; при первом запуске `-m` и `-v` обязательны; `-batch` и `-fs` не запоминаются | `repo_user_area/exe/as_exe.py`, строки 46–68 |
| `README.md`: поддерживаемые площадки (IPP tok, hz-ld-prod, ..., freia) | обычного Linux-ПК среди них нет; `get_platform` возвращает `mypc`, а `platform/env.mypc` отсутствует | `get_platform`; [`DETAILED.md`](../../installer/docs/DETAILED.md), разд. 3.2 |

## 11. Где смотреть в исходниках

| файл | что там |
|---|---|
| `README.md`, `LICENSE` | сопровождающие, лицензия LGPL v2.1+, краткая установка, список площадок |
| `.gitignore` | граница пользовательской области; маркеры слияния `stars3d` |
| `repo_user_area/exe/version` | баннер «Version 8.6, January 2026», авторы |
| `docu/main.tex`, `docu/section12.tex` | основное руководство 2002 г., идея «генератора кодов», история |
| `docu/summary_of_changes_a70.tex`, `docu/addendum_EF.tex` | изменения ASTRA 7.0 |
| `docu/addendum_A8.tex` | изменения ASTRA 8: структура, CLI, синтаксис, `IPEQL`, FEQIS, flight simulator, стелларатор |
| `docu/compilingmain` | как собрать PDF |
| `astra_variables.json` | реестр переменных (14 групп) |
| `src/astra/astra_main.F90` | главная программа: инициализация, сходимость, цикл по времени, выход |
| `src/astra/read_input.f90` | `readInput`, `read_nml` — что и откуда читается во время выполнения |
| `src/astra/json_vars.f90`, `src/astra/json_rw.f90` | чтение реестра, запись/чтение JSON-вывода |
| `pyparse/parser_main.py`, `pyparse/config.py` | вход транслятора, список уравнений и их коэффициентов |
| `postProcess/json2cdf.py` | склейка JSON → CDF |
| `repo_user_area/exe/as_exe.py`, `repo_user_area/exe/Build`, `repo_user_area/exe/Makefile` | запуск и сборка (глава [02](02-build-and-run.md)) |
| `install.sh`, `modules_install.sh`, `modules_inst_loc.sh` | установка ASTRA и внешних модулей |
| `compReg.sh`, `compareRegressions.py`, `regressions/` | регрессионные тесты |
| `memo_stella`, `vmec/` | стеллараторный пример и связь с VMEC |
| `git log` в `astra-src/` | история с 2022-10-13 (`14208069`) |
