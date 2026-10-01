# Документация по исходникам ASTRA 8

Русскоязычная документация по устройству транспортного кода ASTRA 8, собранная по исходникам `astra-src` (ветка `local`, коммит `67281112`, состояние на 2026-10-01) и по оригинальному мануалу `../docu/`.

Вводные для агентов — [`../AGENTS.md`](../AGENTS.md). Установка на Ubuntu/WSL — в монорепо: [`../../installer/README.md`](../../installer/README.md), подробности сборки и патчей — [`../../installer/docs/DETAILED.md`](../../installer/docs/DETAILED.md).

## С чего начать

| если нужно | читать |
|---|---|
| понять, что это за программа и как её используют на токамаке | [10](10-tokamak-usage.md), затем [01](01-overview.md) |
| запустить расчёт и прочитать результат | [11](11-usage-guide.md), [02](02-build-and-run.md), [08](08-graphics-and-output.md) |
| написать или поправить модель | [03](03-model-language.md), [09](09-user-area-catalog.md), [glossary](glossary.md) |
| подготовить данные разряда | [04](04-experiment-data.md), [11](11-usage-guide.md) |
| разобраться в физике и численных схемах | [05](05-solver-core.md), [06](06-equilibrium.md), [07](07-heating-and-external-modules.md) |

## Главы

| № | документ | о чём |
|---|---|---|
| 01 | [`01-overview.md`](01-overview.md) | что такое ASTRA, история A7 → A8, авторы и лицензия, карта всего дерева исходников, путь данных от модели до результата |
| 02 | [`02-build-and-run.md`](02-build-and-run.md) | конвейер сборки и запуска: `install.sh`, `get_platform`, `exe/as_exe(.py)` со всеми опциями, `Build`, Makefile-ы, `astra_rc`, внешние модули `$ASTRA_EXT`, каталоги сборки |
| 03 | [`03-model-language.md`](03-model-language.md) | язык файла модели `equ/<m>` (синтаксис A8, отличия от A7) и Python-парсер `pyparse/`: что и куда он генерирует |
| 04 | [`04-experiment-data.md`](04-experiment-data.md) | файлы эксперимента `exp/<e>`, U-файлы `udb/`, namelist-ы `exp/nml`, `exp/cnf`, `exp/nbi`, путь данных до переменных Fortran |
| 05 | [`05-solver-core.md`](05-solver-core.md) | Fortran-ядро `src/astra`, `src/misc`: главный цикл, уравнения переноса, сетка, численная схема, шаг по времени, IMAS, сообщения об ошибках |
| 06 | [`06-equilibrium.md`](06-equilibrium.md) | равновесие: `IPEQL`, встроенные решатели, FEQIS (свободная граница, катушки, green-матрицы), VMEC, стеллараторный режим |
| 07 | [`07-heating-and-external-modules.md`](07-heating-and-external-modules.md) | нагрев и токовая генерация (NBI, RABBIT, TORBEAM, TORIC, ...), модели турбулентного переноса через IPC (TGLF, NEO, QuaLiKiZ, QLKNN), что доступно на обычной Ubuntu |
| 08 | [`08-graphics-and-output.md`](08-graphics-and-output.md) | X11-интерфейс `src/graph`, вывод json/CDF, `postProcess/`, регрессионные эталоны и их сравнение |
| 09 | [`09-user-area-catalog.md`](09-user-area-catalog.md) | каталог `repo_user_area/`: все модели `equ/`, эксперименты `exp/`, `udb/`, формулы `fml/`, функции `fnc/`, подпрограммы `sbr/` |
| 10 | [`10-tokamak-usage.md`](10-tokamak-usage.md) | как ASTRA используется на токамаке и в лаборатории «Глобус-М2»: интерпретирующий и предсказательный расчёт, входы из эксперимента, что получают физики, литература |
| 11 | [`11-usage-guide.md`](11-usage-guide.md) | практическое руководство: запуск, пошаговый пример, чтение CDF в Python, подготовка разряда, перенос из A7, типичные проблемы |
| — | [`glossary.md`](glossary.md) | переменные ядра, единицы, `astra_variables.json`, термины и сокращения |

## Соглашения

- Пути к файлам даны относительно корня `astra-src`, если не сказано иное. Пути вида `src/tmp/*`, `bin/`, `ncdf_out/` появляются только после сборки (в рабочей установке, например `~/astra/a8`).
- Что не удалось подтвердить по коду, помечено «не проверено» или «предположительно».
- Рабочая установка для сверки — `~/astra/a8` (коммит `e66f6ba1`, отличается от `67281112` правками в `src/misc/imas_ids.f90` и `repo_user_area/exe/astra_rc`).
- Оригинальный мануал (`docu/main.tex`, IPP 5/98, 2002) описывает ASTRA 6/7; изменения A8 — `docu/addendum_A8.tex`. Где мануал расходится с кодом, главы опираются на код и отмечают расхождение.
