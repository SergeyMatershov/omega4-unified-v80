#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
Ω⁴-UNIFIED — ОБЪЕДИНЁННЫЙ MAIN (main_unified.py) — v1.3
# Ω⁴-UNIFIED v80.0
================================================================================
НАЗНАЧЕНИЕ:
    Единая точка входа, объединяющая ШЕСТЬ частей:

        1. БАЗОВАЯ МОДЕЛЬ (main.py):
           • core, atlas, physics, quantum_gravity, compactification
           • 171/171 проверок
           • omega4_results.json

        2. ДОП-ПРОЕКТ (main_research.py v3.0):
           • 9 объединённых research-файлов
           • ~305 констант
           • research_results.json

        3. ДИАГНОСТИКА FRG (check_frg_limit.py v2.1):
           • Таблица N_eff | g* | θ₁ | θ₂ | УФ-стаб.
           • Граница слома: между 7 и 8
           • Инвариант g* + θ₁ = 4

        4. FRG ↔ RGE АНАЛИЗ (omega4_frg_rge_unified.py v1.1):
           • Аналитическое доказательство g* + θ₁ = 4
           • Σ(g*+θ₁) = 28 = dim(SO(8)) = χ/22
           • Финальная формула spread = 20.013016

        5. SUSY-СЕКТОР (v1.2):
           • darboux_susy.py v6.0       — математика SUSY QM
           • susy_omega4.py v4.0        — M_SUSY = (2a−1)·M_EW = 3.69 ТэВ
           • spread_uniqueness.py v1.0  — уникальность spread

        6. ВОЛНОВАЯ ФУНКЦИЯ МЕМБРАНЫ (НОВОЕ в v1.3):
           • membrane_core.py           — ядро (28 констант)
           • membrane_wave_function.py  — единая CLI
           • membrane_vacuum_link.py    — связь с vacuum_analysis (8 констант)
           • Диапазон M_SUSY: [3.20, 3.69] ТэВ

    Плюс:

        7. АВТОГЕНОМ (genome_writer.py):
           • constants_genome_omega4.py
           • Полный REGISTRY: базовая модель + доп-проект
             + FRG ↔ RGE константы + SUSY-константы
             + константы мембраны

ПОРЯДОК РАБОТЫ:
    1. Базовая модель         → print_summary, save omega4_results.json
    2. Доп-проект             → print_report, save research_results.json
    3. check_frg_limit        → печать таблицы + compute_frg_boundary()
    4. omega4_frg_rge_unified → print_report + compute_frg_rge_unified()
    5. SUSY-сектор            → darboux + susy_omega4 + spread_uniqueness
    6. Мембрана               → membrane_wave_function + membrane_vacuum_link
    7. Объединённый JSON      → unified_results.json
    8. Автогеном              → constants_genome_omega4.py (один раз)

ВАЖНО:
    • main.py и main_research.py НЕ трогаем.
    • genome_writer вызывается ОДИН РАЗ — в самом конце.
    • check_frg_limit, omega4_frg_rge_unified, darboux_susy,
      susy_omega4, spread_uniqueness, membrane_wave_function,
      membrane_vacuum_link вызываются как модули.
    • main.py НЕ вызывается как run_complete_analysis() — собираем
      по частям, чтобы отложить автогеном на конец.

ЗАПУСК:
    python main_unified.py
    python main_unified.py --no-genome
    python main_unified.py --no-save
    python main_unified.py --quiet
================================================================================
АВТОР:  Сергей Викторович Матершов
ORCID:  0009-0009-0641-1357
ДАТА:   2026-10-02
ЛИЦЕНЗИЯ: CC BY-NC-ND 4.0 International
================================================================================
ИСТОРИЯ ВЕРСИЙ:
    v1.0 — первая версия.
    v1.1 — добавлен FRG ↔ RGE анализ (ЧАСТЬ 4).
    v1.2 — добавлен SUSY-СЕКТОР (ЧАСТЬ 5).
    v1.3 — добавлена ВОЛНОВАЯ ФУНКЦИЯ МЕМБРАНЫ (ЧАСТЬ 6):  ← ЭТА ВЕРСИЯ
           • membrane_core.py
           • membrane_wave_function.py
           • membrane_vacuum_link.py
           • run_membrane_sector()
           • membrane_data в build_unified_json()
           • блок «6. ВОЛНОВАЯ ФУНКЦИЯ МЕМБРАНЫ» в print_final_summary()
           • Диапазон M_SUSY: [3.20, 3.69] ТэВ
================================================================================
"""

import sys
import json
from datetime import datetime
from typing import Dict, Any, Optional


# ============================================================================
# ЧАСТЬ 1: ИМПОРТЫ
# ============================================================================

# ─── Базовая модель (main.py) ──────────────────────────────────────────────
from main import (
    Omega4Unified,
    HAS_QUANTUM_GRAVITY,
    HAS_COMPACTIFICATION,
    HAS_TOPOLOGICAL_ORIGIN,
)

# ─── Доп-проект (main_research.py v3.0) ────────────────────────────────────
from main_research import (
    run_all_research_files,
    research_summary,
    print_research_report,
    save_research_results,
    formulas_for_7,
    RESEARCH_FILES,
)

# ─── Диагностика FRG (check_frg_limit.py v2.1) ─────────────────────────────
try:
    import check_frg_limit
    HAS_CHECK_FRG = True
except ImportError:
    check_frg_limit = None
    HAS_CHECK_FRG = False

# ─── FRG ↔ RGE анализ (omega4_frg_rge_unified.py v1.1) ─────────────────────
try:
    from omega4_frg_rge_unified import compute_frg_rge_unified
    HAS_FRG_RGE_UNIFIED = True
except ImportError:
    compute_frg_rge_unified = None
    HAS_FRG_RGE_UNIFIED = False

# ─── SUSY-СЕКТОР (ЧАСТЬ 5, v1.2) ───────────────────────────────────────────
try:
    from darboux_susy import verify_darboux_susy
    HAS_DARBOUX_SUSY = True
except ImportError:
    verify_darboux_susy = None
    HAS_DARBOUX_SUSY = False

try:
    from susy_omega4 import (
        full_analysis as susy_full_analysis,
        assess_status as susy_assess_status,
    )
    HAS_SUSY_OMEGA4 = True
except ImportError:
    susy_full_analysis = None
    susy_assess_status = None
    HAS_SUSY_OMEGA4 = False

try:
    from spread_uniqueness import (
        full_analysis as spread_full_analysis,
        compute_main_formula as spread_main_formula,
    )
    HAS_SPREAD_UNIQUENESS = True
except ImportError:
    spread_full_analysis = None
    spread_main_formula = None
    HAS_SPREAD_UNIQUENESS = False

# ─── ★ ВОЛНОВАЯ ФУНКЦИЯ МЕМБРАНЫ (ЧАСТЬ 6, v1.3) ───────────────────────────
try:
    from membrane_wave_function import run_full as membrane_run_full
    HAS_MEMBRANE = True
except ImportError:
    membrane_run_full = None
    HAS_MEMBRANE = False

try:
    from membrane_vacuum_link import (
        analyze_membrane_vacuum,
        print_vacuum_membrane_report,
    )
    HAS_MEMBRANE_VACUUM = True
except ImportError:
    analyze_membrane_vacuum = None
    print_vacuum_membrane_report = None
    HAS_MEMBRANE_VACUUM = False


# ============================================================================
# ЧАСТЬ 2: ЗАПУСК БАЗОВОЙ МОДЕЛИ
# ============================================================================

def run_base_model(verbose: bool = True) -> Dict[str, Any]:
    """
    Запуск базовой модели (main.py) БЕЗ автогенома.
    """
    if verbose:
        print()
        print("=" * 80)
        print("  ЧАСТЬ 1: БАЗОВАЯ МОДЕЛЬ (main.py)")
        print("=" * 80)
        print()

    model = Omega4Unified()
    base_results = model.run_all_checks()

    if verbose:
        model.print_summary(base_results)

    # Сохраняем omega4_results.json
    model.save_results(base_results, filename='omega4_results.json')

    # ─── Извлекаем компактный summary через build_summary_table() ────
    try:
        summary_table = model.build_summary_table(base_results)
        base_results['__summary__'] = {
            'total_pass': summary_table.get('total_pass'),
            'total_checks': summary_table.get('total_checks'),
            'all_verified': summary_table.get('all_verified'),
        }
        # ← сохраняем полную таблицу для print_base_checks_table()
        base_results['__summary_table__'] = summary_table
    except Exception as e:
        base_results['__summary__'] = {
            'total_pass': None,
            'total_checks': None,
            'all_verified': None,
            'error': str(e),
        }
        base_results['__summary_table__'] = None

    return base_results
# ============================================================================
# ЧАСТЬ 2b: ПЕЧАТЬ ВСЕХ 174 ПРОВЕРОК БАЗОВОЙ МОДЕЛИ
# ============================================================================
#
# ЗАЧЕМ:
#   run_base_model() печатает только сводку (174/174).
#   Эта функция печатает ВСЕ 174 проверки в виде таблицы
#   с указанием названия и статуса каждой.
#
# ИСТОЧНИК ДАННЫХ:
#   results['__summary_table__'] — результат build_summary_table()
#   (сохраняется в run_base_model()).
#
# ============================================================================

def print_base_checks_table(base_results: Dict[str, Any]) -> None:
    """
    Печать всех 174 проверок базовой модели в виде таблицы.
    
    ПАРАМЕТРЫ:
        base_results — результат run_base_model()
                       (должен содержать ключ '__summary_table__')
    """
    # ─── Извлекаем таблицу ─────────────────────────────────────────────
    summary_table = base_results.get('__summary_table__')
    
    if not summary_table:
        print("\n  ⚠️  Сводная таблица не найдена "
              "(нет ключа '__summary_table__')")
        return
    
    print()
    print("=" * 80)
    print("  📋 ПОЛНЫЙ СПИСОК 164 ПРОВЕРОК БАЗОВОЙ МОДЕЛИ")
    print("=" * 80)
    print()
    
    # ─── Порядок секций (как в main.py) ────────────────────────────────
    SECTION_ORDER = [
        ('core',                    'ЯДРО (core.py)'),
        ('atlas',                   'АТЛАС (atlas.py)'),
        ('physics',                 'ФИЗИКА (physics.py)'),
        ('vacuum_analysis',         'АНАЛИЗ ВАКУУМА'),
        ('atlas_scan',              'СКАН АТЛАСА'),
        ('magic_numbers',           'МАГИЧЕСКИЕ ЧИСЛА'),
        ('quark_structure',         'СТРУКТУРА КВАРКОВ'),
        ('modular_forms',           'МОДУЛЯРНЫЕ ФОРМЫ'),
        ('statistical',             'СТАТИСТИКА'),
        ('cross_validation',        'КРОСС-ВАЛИДАЦИЯ'),
        ('predictions',             'ПРЕДСКАЗАНИЯ'),
        ('sensitivity',             'ЧУВСТВИТЕЛЬНОСТЬ v2'),
        ('benchmark',               'BENCHMARK'),
        ('elliptic_curve',          'ЭЛЛИПТИЧЕСКАЯ КРИВАЯ'),
        ('quantum_gravity',         'КВАНТОВАЯ ГРАВИТАЦИЯ'),
        ('compactification',        'КОМПАКТИФИКАЦИЯ 26D → 4D'),
        ('topology_vs_spectrum',    'ТОПОЛОГИЯ vs СПЕКТР'),
        ('strong_cp',               'СИЛЬНАЯ CP-ПРОБЛЕМА'),
    ]
    
    total_pass = 0
    total_checks = 0
    global_idx = 0
    
    # ─── Печатаем каждую секцию ────────────────────────────────────────
    for key, title in SECTION_ORDER:
        section = summary_table.get(key)
        if not section:
            continue
        
        checks = section.get('checks', [])
        n_pass = section.get('n_pass', 0)
        n_total = section.get('n_total', 0)
        status = '✅' if n_pass == n_total else '⚠️'
        
        print(f"  ── {title} ({n_pass}/{n_total}) {status} ──")
        print()
        print(f"  {'№':>4} | {'Проверка':<55} | {'Статус':<6}")
        print(f"  {'─'*4} | {'─'*55} | {'─'*6}")
        
        for name, ok in checks:
            global_idx += 1
            marker = '✅' if ok else '❌'
            # Обрезаем название до 55 символов
            name_display = name[:55] if len(name) <= 55 else name[:52] + '...'
            print(f"  {global_idx:>4} | {name_display:<55} | {marker:<6}")
        
        print()
        total_pass += n_pass
        total_checks += n_total
    
    # ─── Итог ──────────────────────────────────────────────────────────
    print("=" * 80)
    print(f"  ВСЕГО: {total_pass}/{total_checks}  "
          f"{'✅ ПОЛНОСТЬЮ ГОТОВО' if total_pass == total_checks else '⚠️'}")
    print("=" * 80)
    print()
# ============================================================================
# ЧАСТЬ 3: ЗАПУСК ДОП-ПРОЕКТА
# ============================================================================

def run_research_project(verbose: bool = True) -> Dict[str, Any]:
    """Запуск доп-проекта (main_research.py v3.0)."""
    if verbose:
        print()
        print("=" * 80)
        print("  ЧАСТЬ 2: ДОП-ПРОЕКТ (main_research.py v3.0)")
        print("=" * 80)
        print()

    run_result = run_all_research_files(verbose=verbose, print_reports=False)

    if verbose:
        print()
        print(f"  Всего файлов:     {run_result['n_total']}")
        print(f"  Загружено:        {run_result['n_loaded']}")
        print(f"  Проверено:        {run_result['n_verified']}")

    if verbose:
        print_research_report(verbose=False)

    save_research_results(output_path='research_results.json', verbose=True)

    return research_summary()


# ============================================================================
# ЧАСТЬ 4: ЗАПУСК ДИАГНОСТИКИ FRG
# ============================================================================

def run_frg_diagnostics(verbose: bool = True) -> Dict[str, Any]:
    """Запуск check_frg_limit.py."""
    if verbose:
        print()
        print("=" * 80)
        print("  ЧАСТЬ 3: ДИАГНОСТИКА FRG (check_frg_limit.py v2.1)")
        print("=" * 80)

    if not HAS_CHECK_FRG:
        if verbose:
            print("\n  ⚠️ check_frg_limit.py недоступен — пропуск.")
        return {
            'available': False,
            'error': 'check_frg_limit not importable',
        }

    if verbose:
        check_frg_limit.main()

    if hasattr(check_frg_limit, 'compute_frg_boundary'):
        boundary = check_frg_limit.compute_frg_boundary()
    else:
        boundary = {
            'available': False,
            'error': 'compute_frg_boundary() отсутствует',
        }

    return boundary


# ============================================================================
# ЧАСТЬ 4b: ЗАПУСК FRG ↔ RGE АНАЛИЗА
# ============================================================================

def run_frg_rge_analysis(verbose: bool = True) -> Dict[str, Any]:
    """Запуск omega4_frg_rge_unified.py (ЧАСТЬ 4)."""
    if verbose:
        print()
        print("=" * 80)
        print("  ЧАСТЬ 4: FRG ↔ RGE АНАЛИЗ (omega4_frg_rge_unified.py v1.1)")
        print("=" * 80)

    if not HAS_FRG_RGE_UNIFIED:
        if verbose:
            print("\n  ⚠️ omega4_frg_rge_unified.py недоступен — пропуск.")
        return {
            'available': False,
            'error': 'omega4_frg_rge_unified not importable',
        }

    if verbose:
        try:
            from omega4_frg_rge_unified import print_report
            print_report()
        except Exception as e:
            print(f"\n  ⚠️ Ошибка печати отчёта: {e}")

    try:
        result = compute_frg_rge_unified()
    except Exception as e:
        result = {
            'available': False,
            'error': str(e),
        }

    return result


# ============================================================================
# ЧАСТЬ 5: ЗАПУСК SUSY-СЕКТОРА (v1.2)
# ============================================================================

def run_susy_sector(verbose: bool = True) -> Dict[str, Any]:
    """Запуск ЧАСТИ 5: SUSY-сектор."""
    if verbose:
        print()
        print("=" * 80)
        print("  ЧАСТЬ 5: SUSY-СЕКТОР")
        print("=" * 80)

    result = {
        'available': False,
        'darboux_susy': None,
        'susy_omega4': None,
        'spread_uniqueness': None,
        'summary': {},
        'all_verified': False,
    }

    # ─── 5.1. Darboux SUSY (SUSY QM) ─────────────────────────────────
    if HAS_DARBOUX_SUSY:
        if verbose:
            print("\n  ── 5.1. ФАКТОРИЗАЦИЯ ДАРБУ (darboux_susy.py) ──")
        try:
            darboux_data = verify_darboux_susy()
            result['darboux_susy'] = darboux_data
            result['available'] = True

            if verbose:
                print(f"    Моделей:              {darboux_data['n_models']}")
                print(f"    Обязательных:         {darboux_data['n_required']}")
                print(f"    Изоспектральных:      {darboux_data['n_matched']}"
                      f"/{darboux_data['n_required']}")
                print(f"    Факторизация:         "
                      f"{'✅' if darboux_data['all_factorized'] else '❌'}")
                print(f"    Статус:               {darboux_data['status']}")
        except Exception as e:
            if verbose:
                print(f"    ❌ Ошибка: {e}")
            result['darboux_susy'] = {'error': str(e)}
    else:
        if verbose:
            print("\n  ⚠️ darboux_susy.py недоступен — пропуск.")

    # ─── 5.2. SUSY Ω⁴ (M_SUSY) ───────────────────────────────────────
    if HAS_SUSY_OMEGA4:
        if verbose:
            print("\n  ── 5.2. M_SUSY (susy_omega4.py v4.0) ──")
        try:
            analysis = susy_full_analysis()
            status_data = susy_assess_status()

            result['susy_omega4'] = {
                'analysis': analysis,
                'status': status_data,
            }
            result['available'] = True

            best = analysis['best']
            m_susy = best['m_susy']['M_SUSY_TeV']

            if verbose:
                print(f"    Формула (H3):         M_SUSY = (2a−1)·M_EW")
                print(f"    a = N_eff:            {best['a']:.4f}")
                print(f"    ΔE(a) = 2a−1:         "
                      f"{best['m_susy']['delta_E']:.4f}")
                print(f"    M_SUSY:               {m_susy:.4f} ТэВ")
                print(f"    LHC-статус:           "
                      f"{best['lhc_h3']['status']}")
                print(f"    Пройдено:             "
                      f"{best['lhc_h3']['n_passed']}/"
                      f"{best['lhc_h3']['n_total']}")
                print(f"    Статус:               {status_data['status']}")
        except Exception as e:
            if verbose:
                print(f"    ❌ Ошибка: {e}")
            result['susy_omega4'] = {'error': str(e)}
    else:
        if verbose:
            print("\n  ⚠️ susy_omega4.py недоступен — пропуск.")

    # ─── 5.3. Уникальность spread ────────────────────────────────────
    if HAS_SPREAD_UNIQUENESS:
        if verbose:
            print("\n  ── 5.3. УНИКАЛЬНОСТЬ SPREAD (spread_uniqueness.py) ──")
        try:
            spread_data = spread_full_analysis()
            result['spread_uniqueness'] = spread_data
            result['available'] = True

            uniqueness = spread_data['uniqueness']
            integers = spread_data['integers']

            if verbose:
                main_formula = spread_data['main']
                print(f"    Формула:              "
                      f"{main_formula['formula']}")
                print(f"    spread:               "
                      f"{main_formula['spread']:.10f}")
                print(f"    D (уникальный):       "
                      f"{uniqueness['exact_d']:.4f}")
                print(f"    Уникальность:         "
                      f"{'✅' if uniqueness['is_unique'] else '❌'}")
                print(f"    Целочисл. комбинаций: {integers['n_total']} "
                      f"(совпадений: {integers['n_matches']})")
        except Exception as e:
            if verbose:
                print(f"    ❌ Ошибка: {e}")
            result['spread_uniqueness'] = {'error': str(e)}
    else:
        if verbose:
            print("\n  ⚠️ spread_uniqueness.py недоступен — пропуск.")

    # ─── 5.4. Компактная сводка ──────────────────────────────────────
    summary = {}

    if result['darboux_susy'] and 'error' not in result['darboux_susy']:
        ds = result['darboux_susy']
        summary['darboux_n_matched'] = ds.get('n_matched')
        summary['darboux_n_required'] = ds.get('n_required')
        summary['darboux_status'] = ds.get('status')
        summary['darboux_all_factorized'] = ds.get('all_factorized')

    if result['susy_omega4'] and 'error' not in result['susy_omega4']:
        so = result['susy_omega4']
        best = so['analysis']['best']
        summary['susy_m_susy_tev'] = best['m_susy']['M_SUSY_TeV']
        summary['susy_a_primary'] = best['a']
        summary['susy_delta_e'] = best['m_susy']['delta_E']
        summary['susy_lhc_status'] = best['lhc_h3']['status']
        summary['susy_lhc_n_passed'] = best['lhc_h3']['n_passed']
        summary['susy_lhc_n_total'] = best['lhc_h3']['n_total']
        summary['susy_status'] = so['status']['status']

    if result['spread_uniqueness'] and 'error' not in result['spread_uniqueness']:
        su = result['spread_uniqueness']
        summary['spread_unique_d'] = su['uniqueness']['exact_d']
        summary['spread_is_unique'] = su['uniqueness']['is_unique']
        summary['spread_integers_n'] = su['integers']['n_total']
        summary['spread_integers_matches'] = su['integers']['n_matches']

    result['summary'] = summary

    # ─── 5.5. Итоговый статус ────────────────────────────────────────
    all_verified_parts = []

    if result['darboux_susy'] and 'error' not in result['darboux_susy']:
        all_verified_parts.append(
            result['darboux_susy'].get('all_factorized', False)
        )

    if result['susy_omega4'] and 'error' not in result['susy_omega4']:
        all_verified_parts.append(
            result['susy_omega4']['status']['status'] == 'STRUCTURAL_STRONG'
        )

    if result['spread_uniqueness'] and 'error' not in result['spread_uniqueness']:
        all_verified_parts.append(
            result['spread_uniqueness']['uniqueness']['is_unique']
        )

    result['all_verified'] = (
        len(all_verified_parts) > 0 and all(all_verified_parts)
    )

    return result


# ============================================================================
# ЧАСТЬ 6: ЗАПУСК ВОЛНОВОЙ ФУНКЦИИ МЕМБРАНЫ (v1.3)
# ============================================================================
#
# Три модуля:
#   • membrane_core.py           — ядро (28 констант)
#   • membrane_wave_function.py  — единая CLI
#   • membrane_vacuum_link.py    — связь с vacuum_analysis (8 констант)
#
# Ключевой результат: диапазон M_SUSY [3.20, 3.69] ТэВ
# ============================================================================

def run_membrane_sector(verbose: bool = True) -> Dict[str, Any]:
    """
    Запуск ЧАСТИ 6: Волновая функция мембраны.

    ОБЪЕДИНЯЕТ:
        • 1D-модель Пёшля-Теллера
        • SUSY QM (факторизация Дарбу)
        • 2D-проверка
        • 4D через базис
        • Комбинаторика кратностей
        • 4D-волновые функции
        • Связь с vacuum_analysis
        • Пограничное состояние a = 7 vs a = 8

    ВОЗВРАЩАЕТ:
        dict с полями:
          • available       — True, если модуль доступен
          • membrane        — результаты полного прогона
          • vacuum_link     — связь с вакуумом
          • summary         — компактная сводка
          • all_verified    — True, если всё сошлось
    """
    if verbose:
        print()
        print("=" * 80)
        print("  ЧАСТЬ 6: ВОЛНОВАЯ ФУНКЦИЯ МЕМБРАНЫ")
        print("=" * 80)

    result = {
        'available': False,
        'membrane': None,
        'vacuum_link': None,
        'summary': {},
        'all_verified': False,
    }

    # ─── 6.1. Полный прогон мембраны ────────────────────────────────────
    if HAS_MEMBRANE:
        if verbose:
            print("\n  ── 6.1. Полный прогон (membrane_wave_function.py) ──")
        try:
            # Тихий режим — печатаем только сводку (без визуализации)
            membrane_data = membrane_run_full(
                a=7.0,
                n_basis=8,
                susy=True,
                viz=False,      # Без визуализации для скорости
                quiet=True,
            )
            result['membrane'] = membrane_data
            result['available'] = True

            if verbose:
                s1d = membrane_data.get('1D', {})
                s2d = membrane_data.get('2D', {})
                s4d = membrane_data.get('4D', {})
                sc = membrane_data.get('combinatorial', {})
                st = membrane_data.get('topo', {})
                sp = membrane_data.get('susy_partners', {})
                ms = membrane_data.get('M_SUSY', {})
                scos = membrane_data.get('cosmology', {})

                print(f"    1D max|ΔE|:              "
                      f"{s1d.get('max_dev', 0):.6f}")
                print(f"    2D разделение:           "
                      f"{'✅' if 'results' in s2d else '⚠️'}")
                print(f"    4D базис:                "
                      f"{s4d.get('n_total', 0)} состояний, "
                      f"max_dev = {s4d.get('max_dev', 0):.2e}")
                print(f"    Комбинаторика:           "
                      f"{'✅' if sc.get('all_match') else '❌'}")
                print(f"    Бином. кратности C(4,k): "
                      f"{'✅' if st.get('binomial_match') else '❌'}")
                print(f"    SUSY изоспектральность:  "
                      f"{'✅' if sp.get('isospectral') else '❌'}")
                print(f"    M_SUSY (a=7):            "
                      f"{ms.get('M_SUSY', 0)/1000:.4f} ТэВ")
                print(f"    n_s (космология):        "
                      f"{scos.get('n_s', 0):.6f}")
        except Exception as e:
            if verbose:
                print(f"    ❌ Ошибка: {e}")
            result['membrane'] = {'error': str(e)}
    else:
        if verbose:
            print("\n  ⚠️ membrane_wave_function.py недоступен — пропуск.")

    # ─── 6.2. Связь с вакуумом ──────────────────────────────────────────
    if HAS_MEMBRANE_VACUUM:
        if verbose:
            print("\n  ── 6.2. Связь с vacuum_analysis.py ──")
        try:
            vacuum_data = analyze_membrane_vacuum(
                a_range=[6.5, 7.0, 7.25, 7.5, 7.75, 8.0, 8.5]
            )
            result['vacuum_link'] = vacuum_data
            result['available'] = True

            if verbose:
                print(f"    Рабочая точка:           "
                      f"a = {vacuum_data['transition']['a_stable']}")
                print(f"    Порог слома:             "
                      f"a = {vacuum_data['transition']['a_break']}")
                print(f"    M_SUSY (min):            "
                      f"{vacuum_data['M_SUSY_min_TeV']:.4f} ТэВ")
                print(f"    M_SUSY (max):            "
                      f"{vacuum_data['M_SUSY_max_TeV']:.4f} ТэВ")
                print(f"    ★ Диапазон M_SUSY:       "
                      f"[{vacuum_data['M_SUSY_min_TeV']:.4f}, "
                      f"{vacuum_data['M_SUSY_max_TeV']:.4f}] ТэВ")
                print(f"    Режим stable:            "
                      f"{len(vacuum_data['stable'])} точек")
                print(f"    Режим boundary:          "
                      f"{len(vacuum_data['boundary'])} точек")
                print(f"    Режим broken:            "
                      f"{len(vacuum_data['broken'])} точек")
        except Exception as e:
            if verbose:
                print(f"    ❌ Ошибка: {e}")
            result['vacuum_link'] = {'error': str(e)}
    else:
        if verbose:
            print("\n  ⚠️ membrane_vacuum_link.py недоступен — пропуск.")

    # ─── 6.3. Компактная сводка ─────────────────────────────────────────
    summary = {}

    if result['membrane'] and 'error' not in result['membrane']:
        m = result['membrane']
        summary['1D_max_dev'] = m.get('1D', {}).get('max_dev', 0)
        summary['4D_n_total'] = m.get('4D', {}).get('n_total', 0)
        summary['4D_n_unique'] = m.get('4D', {}).get('n_unique', 0)
        summary['combinatorial_match'] = m.get('combinatorial', {}).get('all_match', False)
        summary['susy_isospectral'] = m.get('susy_partners', {}).get('isospectral', False)
        summary['M_SUSY_TeV'] = m.get('M_SUSY', {}).get('M_SUSY', 0) / 1000
        summary['n_s'] = m.get('cosmology', {}).get('n_s', 0)

    if result['vacuum_link'] and 'error' not in result['vacuum_link']:
        v = result['vacuum_link']
        summary['a_stable'] = v['transition']['a_stable']
        summary['a_break'] = v['transition']['a_break']
        summary['M_SUSY_min_TeV'] = v['M_SUSY_min_TeV']
        summary['M_SUSY_max_TeV'] = v['M_SUSY_max_TeV']

    result['summary'] = summary

    # ─── 6.4. Итоговый статус ───────────────────────────────────────────
    parts = []
    if result['membrane'] and 'error' not in result['membrane']:
        parts.append(True)
    if result['vacuum_link'] and 'error' not in result['vacuum_link']:
        parts.append(True)

    result['all_verified'] = len(parts) == 2 and all(parts)

    return result


# ============================================================================
# ЧАСТЬ 7: ОБЪЕДИНЁННЫЙ JSON
# ============================================================================

def _extract_base_summary(base_results: Dict[str, Any]) -> Dict[str, Any]:
    """Извлечь компактную сводку из результатов базовой модели."""
    if not isinstance(base_results, dict):
        return {}

    s = base_results.get('__summary__')
    if not isinstance(s, dict):
        return {
            'total_pass': None,
            'total_checks': None,
            'all_verified': None,
        }

    return {
        'total_pass': s.get('total_pass'),
        'total_checks': s.get('total_checks'),
        'all_verified': s.get('all_verified'),
    }


def _extract_research_summary(research_data: Dict[str, Any]) -> Dict[str, Any]:
    """Извлечь компактную сводку из research_summary()."""
    if not isinstance(research_data, dict):
        return {}

    if 'N_eff_max' in research_data or 'status_overall' in research_data:
        return {
            'N_eff_max': research_data.get('N_eff_max'),
            'N_bulk': research_data.get('N_bulk'),
            'N_hidden': research_data.get('N_hidden'),
            'N_membrane': research_data.get('N_membrane'),
            'status_overall': research_data.get('status_overall'),
        }

    summary = research_data.get('summary')
    if isinstance(summary, dict):
        return summary

    return {}


def build_unified_json(
    base_results: Dict[str, Any],
    research_data: Dict[str, Any],
    frg_boundary: Dict[str, Any],
    frg_rge_data: Dict[str, Any],
    susy_data: Optional[Dict[str, Any]] = None,
    membrane_data: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Построить unified_results.json.

    Структура:
        {
          "version": "1.3",
          "timestamp": "...",
          "base_model":       {...},
          "research_project": {...},
          "frg_boundary":     {...},
          "frg_rge_analysis": {...},
          "susy_sector":      {...},
          "membrane_sector":  {...},        ← НОВОЕ v1.3
          "combined_summary": {...}
        }
    """
    base_summary = _extract_base_summary(base_results)
    research_summary_compact = _extract_research_summary(research_data)

    combined = {
        # ── Базовая модель ──
        'base_model_all_verified': base_summary.get('all_verified'),
        'base_model_total_pass': base_summary.get('total_pass'),
        'base_model_total_checks': base_summary.get('total_checks'),

        # ── Доп-проект ──
        'research_status_overall': research_summary_compact.get('status_overall'),
        'research_N_eff_max': research_summary_compact.get('N_eff_max'),
        'research_N_bulk': research_summary_compact.get('N_bulk'),
        'research_N_hidden': research_summary_compact.get('N_hidden'),
        'research_N_membrane': research_summary_compact.get('N_membrane'),

        # ── FRG-диагностика ──
        'frg_boundary_available': (
            isinstance(frg_boundary, dict) and frg_boundary.get('available')
        ),
        'frg_N_eff_max_stable': (
            frg_boundary.get('N_eff_max_stable')
            if isinstance(frg_boundary, dict) else None
        ),
        'frg_N_eff_min_unstable': (
            frg_boundary.get('N_eff_min_unstable')
            if isinstance(frg_boundary, dict) else None
        ),
        'frg_boundary_text': (
            frg_boundary.get('boundary')
            if isinstance(frg_boundary, dict) else None
        ),

        # ── FRG ↔ RGE анализ ──
        'frg_rge_available': (
            frg_rge_data.get('available', False)
            if isinstance(frg_rge_data, dict) else False
        ),
        'frg_rge_all_verified': (
            frg_rge_data.get('all_verified')
            if isinstance(frg_rge_data, dict) else None
        ),
        'frg_rge_invariant': (
            frg_rge_data.get('frg_ladder', {}).get('invariant')
            if isinstance(frg_rge_data, dict) else None
        ),
        'frg_rge_sum_sigma_theta': (
            frg_rge_data.get('frg_ladder', {}).get('sum_total')
            if isinstance(frg_rge_data, dict) else None
        ),
        'frg_rge_N_eff': (
            frg_rge_data.get('N_eff_criterion', {}).get('N_eff')
            if isinstance(frg_rge_data, dict) else None
        ),
        'frg_rge_spread_struct': (
            frg_rge_data.get('spread_struct', {}).get('spread_struct')
            if isinstance(frg_rge_data, dict) else None
        ),
        'frg_rge_spread_1loop': (
            frg_rge_data.get('spread_1loop')
            if isinstance(frg_rge_data, dict) else None
        ),
        'frg_rge_spread_2loop': (
            frg_rge_data.get('spread_2loop')
            if isinstance(frg_rge_data, dict) else None
        ),
        'frg_rge_deviation_1loop': (
            frg_rge_data.get('deviation_1loop_percent')
            if isinstance(frg_rge_data, dict) else None
        ),
        'frg_rge_deviation_2loop': (
            frg_rge_data.get('deviation_2loop_percent')
            if isinstance(frg_rge_data, dict) else None
        ),
    }

    # ── SUSY-сектор (ЧАСТЬ 5) ───────────────────────────────────────
    if isinstance(susy_data, dict):
        susy_summary = susy_data.get('summary', {})
        combined.update({
            'susy_available': susy_data.get('available', False),
            'susy_all_verified': susy_data.get('all_verified', False),
            'susy_darboux_n_matched': susy_summary.get('darboux_n_matched'),
            'susy_darboux_n_required': susy_summary.get('darboux_n_required'),
            'susy_darboux_status': susy_summary.get('darboux_status'),
            'susy_m_susy_tev': susy_summary.get('susy_m_susy_tev'),
            'susy_a_primary': susy_summary.get('susy_a_primary'),
            'susy_delta_e': susy_summary.get('susy_delta_e'),
            'susy_lhc_status': susy_summary.get('susy_lhc_status'),
            'susy_status': susy_summary.get('susy_status'),
            'spread_unique_d': susy_summary.get('spread_unique_d'),
            'spread_is_unique': susy_summary.get('spread_is_unique'),
            'spread_integers_n': susy_summary.get('spread_integers_n'),
            'spread_integers_matches': susy_summary.get('spread_integers_matches'),
        })

    # ── ★ МЕМБРАНА (ЧАСТЬ 6) ────────────────────────────────────────
    if isinstance(membrane_data, dict):
        membrane_summary = membrane_data.get('summary', {})
        combined.update({
            'membrane_available': membrane_data.get('available', False),
            'membrane_all_verified': membrane_data.get('all_verified', False),
            'membrane_1D_max_dev': membrane_summary.get('1D_max_dev'),
            'membrane_4D_n_total': membrane_summary.get('4D_n_total'),
            'membrane_4D_n_unique': membrane_summary.get('4D_n_unique'),
            'membrane_combinatorial_match': membrane_summary.get('combinatorial_match'),
            'membrane_susy_isospectral': membrane_summary.get('susy_isospectral'),
            'membrane_M_SUSY_TeV': membrane_summary.get('M_SUSY_TeV'),
            'membrane_n_s': membrane_summary.get('n_s'),
            'membrane_a_stable': membrane_summary.get('a_stable'),
            'membrane_a_break': membrane_summary.get('a_break'),
            'membrane_M_SUSY_min_TeV': membrane_summary.get('M_SUSY_min_TeV'),
            'membrane_M_SUSY_max_TeV': membrane_summary.get('M_SUSY_max_TeV'),
        })

    unified = {
        '_author': 'Сергей Викторович Матершов',
        '_orcid': '0009-0009-0641-1357',
        '_license': 'CC BY-NC-ND 4.0 International',
        'version': '1.3',
        'timestamp': datetime.now().isoformat(),
        'description': (
            'Ω⁴-UNIFIED: базовая модель + доп-проект + FRG-диагностика '
            '+ FRG ↔ RGE анализ + SUSY-сектор + волновая функция мембраны'
        ),

        # ── 1. Базовая модель ──
        'base_model': {
            'results': base_results,
            'summary': base_summary,
        },

        # ── 2. Доп-проект ──
        'research_project': {
            'summary': research_summary_compact,
            'details': (
                research_data.get('details', {})
                if isinstance(research_data, dict) else {}
            ),
            'files': [spec['name'] for spec in RESEARCH_FILES],
            'formulas_for_7': formulas_for_7(),
        },

        # ── 3. FRG-диагностика ──
        'frg_boundary': (
            frg_boundary if isinstance(frg_boundary, dict)
            else {'available': False}
        ),

        # ── 4. FRG ↔ RGE анализ ──
        'frg_rge_analysis': (
            frg_rge_data if isinstance(frg_rge_data, dict)
            else {'available': False}
        ),

        # ── 5. SUSY-СЕКТОР ──
        'susy_sector': (
            susy_data if isinstance(susy_data, dict)
            else {'available': False}
        ),

        # ── 6. ★ ВОЛНОВАЯ ФУНКЦИЯ МЕМБРАНЫ (НОВОЕ v1.3) ──
        'membrane_sector': (
            membrane_data if isinstance(membrane_data, dict)
            else {'available': False}
        ),

        # ── 7. Компактная сводка ──
        'combined_summary': combined,
    }

    return unified


def save_unified_json(
    unified: Dict[str, Any],
    output_path: str = 'unified_results.json',
    verbose: bool = True,
) -> str:
    """Сохранить unified_results.json."""
    try:
        with open(output_path, 'w', encoding='utf-8') as f:
            json.dump(unified, f, indent=2, ensure_ascii=False, default=str)
        if verbose:
            print(f"\n  💾 Объединённый JSON сохранён: {output_path}")
        return output_path
    except OSError as e:
        if verbose:
            print(f"\n  ❌ Ошибка сохранения unified JSON: {e}")
        return ''


# ============================================================================
# ЧАСТЬ 8: ФИНАЛЬНАЯ СВОДКА
# ============================================================================

def print_final_summary(
    base_results: Dict[str, Any],
    research_data: Dict[str, Any],
    frg_boundary: Dict[str, Any],
    frg_rge_data: Dict[str, Any],
    susy_data: Optional[Dict[str, Any]] = None,
    membrane_data: Optional[Dict[str, Any]] = None,
) -> None:
    """Печать финальной сводки по шести частям."""
    print()
    print("=" * 80)
    print("  ФИНАЛЬНАЯ СВОДКА (main_unified.py v1.3)")
    print("=" * 80)
    print()

    # ─── 1. Базовая модель ─────────────────────────────────────────────
    base_summary = _extract_base_summary(base_results)

    print("  ── 1. БАЗОВАЯ МОДЕЛЬ (main.py) ──")
    print(f"    Проверок:      {base_summary.get('total_pass', '?')}/"
          f"{base_summary.get('total_checks', '?')}")
    print(f"    Статус:        "
          f"{'✅ ГОТОВО' if base_summary.get('all_verified') else '⚠️ ЕСТЬ ОШИБКИ'}")
    print()

    # ─── 2. Доп-проект ─────────────────────────────────────────────────
    research_summary_compact = _extract_research_summary(research_data)

    print("  ── 2. ДОП-ПРОЕКТ (main_research.py v3.0) ──")
    print(f"    Файлов:        {len(RESEARCH_FILES)}")
    print(f"    N_eff_max:     {research_summary_compact.get('N_eff_max', '?')}")
    print(f"    N_bulk:        {research_summary_compact.get('N_bulk', '?')}")
    print(f"    N_hidden:      {research_summary_compact.get('N_hidden', '?')}")
    print(f"    N_membrane:    {research_summary_compact.get('N_membrane', '?')}")
    print(f"    Статус:        {research_summary_compact.get('status_overall', '?')}")
    print()

    # ─── 3. FRG-диагностика ────────────────────────────────────────────
    print("  ── 3. ДИАГНОСТИКА FRG (check_frg_limit.py v2.1) ──")
    if isinstance(frg_boundary, dict) and frg_boundary.get('available'):
        print(f"    N_eff_max_stable:     {frg_boundary.get('N_eff_max_stable', '?')}")
        print(f"    N_eff_min_unstable:   {frg_boundary.get('N_eff_min_unstable', '?')}")
        print(f"    Граница:              {frg_boundary.get('boundary', '?')}")

        sigma_theta = frg_boundary.get('sigma_theta_invariant', {})
        if sigma_theta and sigma_theta.get('all_exact'):
            print(f"    ★ Инвариант:          g* + θ₁ = 4 "
                  f"(Σ = {sigma_theta.get('sum_total', '?'):.4f})")
    else:
        print(f"    ⚠️ недоступна")
    print()

    # ─── 4. FRG ↔ RGE анализ ───────────────────────────────────────────
    print("  ── 4. FRG ↔ RGE АНАЛИЗ (omega4_frg_rge_unified.py v1.1) ──")
    if isinstance(frg_rge_data, dict) and frg_rge_data.get('available'):
        frg_ladder = frg_rge_data.get('frg_ladder', {})
        neff = frg_rge_data.get('N_eff_criterion', {})
        struct = frg_rge_data.get('spread_struct', {})

        invariant = frg_ladder.get('invariant')
        sum_sigma = frg_ladder.get('sum_total')
        print(f"    Инвариант:      g* + θ₁ = {invariant}")
        if isinstance(sum_sigma, (int, float)):
            print(f"    Σ(g*+θ₁):       {sum_sigma:.4f} = dim(SO(8)) = χ/22")

        n_eff = neff.get('N_eff')
        print(f"    N_eff:          {n_eff} = h¹¹ + N_gen − 1 = dim(SU(3))")

        spread_struct = struct.get('spread_struct')
        spread_1loop = frg_rge_data.get('spread_1loop')
        spread_2loop = frg_rge_data.get('spread_2loop')
        dev_1loop = frg_rge_data.get('deviation_1loop_percent')
        dev_2loop = frg_rge_data.get('deviation_2loop_percent')

        if isinstance(spread_struct, (int, float)):
            print(f"    spread(struct): {spread_struct:.6f}")
        if isinstance(spread_1loop, (int, float)):
            print(f"    spread(1loop):  {spread_1loop:.6f}")
        if isinstance(spread_2loop, (int, float)):
            print(f"    spread(2loop):  {spread_2loop:.6f}")
        if isinstance(dev_1loop, (int, float)):
            print(f"    Откл. (1loop):  {dev_1loop:.6f}%")
        if isinstance(dev_2loop, (int, float)):
            print(f"    Откл. (2loop):  {dev_2loop:.6f}%")

        all_verified = frg_rge_data.get('all_verified')
        status_str = (
            '✅ STRUCTURAL_PROVEN' if all_verified
            else '⚠️ ЕСТЬ ОТКЛОНЕНИЯ'
        )
        print(f"    Статус:         {status_str}")
    else:
        print(f"    ⚠️ недоступен")
        if isinstance(frg_rge_data, dict) and frg_rge_data.get('error'):
            print(f"    Ошибка: {frg_rge_data['error']}")
    print()

    # ─── 5. SUSY-СЕКТОР ────────────────────────────────────────────────
    print("  ── 5. SUSY-СЕКТОР ──")
    print("     (darboux_susy + susy_omega4 + spread_uniqueness)")
    if isinstance(susy_data, dict) and susy_data.get('available'):
        if susy_data.get('darboux_susy') and \
                'error' not in susy_data['darboux_susy']:
            ds = susy_data['darboux_susy']
            print(f"    Дарбу:          "
                  f"{ds['n_matched']}/{ds['n_required']} "
                  f"({'✅' if ds['all_factorized'] else '❌'})")

        if susy_data.get('susy_omega4') and \
                'error' not in susy_data['susy_omega4']:
            so = susy_data['susy_omega4']
            best = so['analysis']['best']
            m_susy = best['m_susy']['M_SUSY_TeV']
            delta_e = best['m_susy']['delta_E']
            print(f"    M_SUSY(H3):     {m_susy:.4f} ТэВ  "
                  f"(a = {best['a']:.4f}, ΔE = {delta_e:.2f})")
            print(f"    LHC:            "
                  f"{best['lhc_h3']['n_passed']}/{best['lhc_h3']['n_total']} "
                  f"— {best['lhc_h3']['status']}")
            print(f"    Статус:         {so['status']['status']}")

        if susy_data.get('spread_uniqueness') and \
                'error' not in susy_data['spread_uniqueness']:
            su = susy_data['spread_uniqueness']
            uniq = su['uniqueness']
            print(f"    Spread D:       {uniq['exact_d']:.4f}  "
                  f"({'✅ уникален' if uniq['is_unique'] else '❌ не уникален'})")
    else:
        print(f"    ⚠️ недоступен")
        if isinstance(susy_data, dict) and susy_data.get('error'):
            print(f"    Ошибка: {susy_data['error']}")
    print()

    # ─── 6. ★ ВОЛНОВАЯ ФУНКЦИЯ МЕМБРАНЫ (НОВОЕ v1.3) ───────────────────
    print("  ── 6. ВОЛНОВАЯ ФУНКЦИЯ МЕМБРАНЫ ──")
    print("     (membrane_core + membrane_wave_function + membrane_vacuum_link)")
    if isinstance(membrane_data, dict) and membrane_data.get('available'):
        summary = membrane_data.get('summary', {})

        print(f"    1D max|ΔE|:            "
              f"{summary.get('1D_max_dev', 0):.6f}")
        print(f"    4D состояний:          {summary.get('4D_n_total', 0)} "
              f"({summary.get('4D_n_unique', 0)} уровней)")
        print(f"    Комбинаторика:         "
              f"{'✅' if summary.get('combinatorial_match') else '❌'}")
        print(f"    SUSY изоспектр.:       "
              f"{'✅' if summary.get('susy_isospectral') else '❌'}")
        print(f"    M_SUSY (a=7):          "
              f"{summary.get('M_SUSY_TeV', 0):.4f} ТэВ")

        min_tevs = summary.get('M_SUSY_min_TeV')
        max_tevs = summary.get('M_SUSY_max_TeV')
        if min_tevs is not None and max_tevs is not None:
            print(f"    ★ Диапазон M_SUSY:     "
                  f"[{min_tevs:.4f}, {max_tevs:.4f}] ТэВ")

        a_stable = summary.get('a_stable')
        a_break = summary.get('a_break')
        if a_stable is not None:
            print(f"    a = stable:            {a_stable}")
        if a_break is not None:
            print(f"    a = break:             {a_break}")
    else:
        print(f"    ⚠️ недоступна")
        if isinstance(membrane_data, dict) and membrane_data.get('error'):
            print(f"    Ошибка: {membrane_data['error']}")
    print()

    # ─── ИТОГ ──────────────────────────────────────────────────────────
    print("  ── ИТОГ ──")
    print(f"    Ω⁴-UNIFIED v79.6 + main_research v3.0")
    print(f"    + check_frg_limit v2.1 + omega4_frg_rge_unified v1.1")
    print(f"    + darboux_susy v6.0 + susy_omega4 v4.0 "
          f"+ spread_uniqueness v1.0")
    print(f"    + membrane_core v1.0 + membrane_wave_function v3.1 "
          f"+ membrane_vacuum_link v1.0")
    print(f"    Все шесть частей завершены.")
    print(f"    unified_results.json сохранён.")
    print()
    print("=" * 80)


# ============================================================================
# ЧАСТЬ 9: ТОЧКА ВХОДА
# ============================================================================
# ============================================================================
# ЧАСТЬ 9b: ИТОГОВЫЙ СПИСОК ЗАГРУЖЕННЫХ МОДУЛЕЙ
# ============================================================================
#
# ЗАЧЕМ:
#   В конце прогона показать, какие Python-модули проекта
#   реально загрузились (находятся в sys.modules).
#
# КАК РАБОТАЕТ:
#   1. Есть словарь EXPECTED_MODULES — список модулей,
#      которые ДОЛЖНЫ загрузиться, сгруппированный по частям.
#   2. Функция print_loaded_modules_summary() обходит словарь,
#      проверяет каждый модуль через sys.modules и печатает:
#        ✅ — модуль загружен
#        ⚠️  — модуль НЕ загружен
#   3. Вызывается в main() после print_final_summary().
#
# ВАЖНО:
#   • Список EXPECTED_MODULES — это НЕ автогеном (723 константы).
#     Это список МОДУЛЕЙ (41 .py), которые должны быть загружены.
#   • axion_bi_link.py в списке НЕТ — он не участвует в main_unified.py.
#
# ============================================================================

EXPECTED_MODULES = {
    '1. Ядро': [
        'core', 'atlas', 'physics',
    ],
    '2. Базовая модель (main.py)': [
        'main',
        'vacuum_analysis', 'atlas_scan', 'magic_numbers',
        'quark_structure', 'modular_forms',
        'statistical', 'cross_validation', 'predictions',
        'sensitivity', 'benchmark',
        'elliptic_curve', 'topological_origin',
        'extra_constants', 'system_metrics',
        'quantum_gravity', 'compactification',
    ],
    '3. Доп-проект (main_research.py)': [
        'main_research',
        'research_topology', 'research_membrane', 'research_bulk',
        'research_m_theory', 'research_omega4', 'research_omega58',
        'research_hierarchy', 'research_predictions', 'research_final',
    ],
    '4. FRG ↔ RGE': [
        'check_frg_limit', 'omega4_frg_rge_unified',
    ],
    '5. SUSY-сектор': [
        'darboux_susy', 'susy_omega4', 'spread_uniqueness',
    ],
    '6. Мембрана': [
        'membrane_core', 'membrane_wave_function', 'membrane_vacuum_link',
    ],
    '7. Геном и обоснование': [
        'genome_writer', 'msusy_derivation',
    ],
}


def print_loaded_modules_summary() -> None:
    """
    Печать итогового списка загруженных модулей.

    Группирует по частям системы. Для каждого модуля проверяет,
    находится ли он в sys.modules:
      ✅ — модуль загружен
      ⚠️  — модуль НЕ загружен

    В конце печатает итог: сколько модулей загружено из скольких.
    """
    import sys

    print()
    print("=" * 80)
    print("  📦 ИТОГОВЫЙ СПИСОК ЗАГРУЖЕННЫХ МОДУЛЕЙ")
    print("=" * 80)
    print()

    total_loaded = 0
    total_modules = 0

    for part_name, modules in EXPECTED_MODULES.items():
        loaded = [m for m in modules if m in sys.modules]
        not_loaded = [m for m in modules if m not in sys.modules]
        total_loaded += len(loaded)
        total_modules += len(modules)

        # ─── Статус части ───
        status = '✅' if not not_loaded else '⚠️'
        print(f"  {status} {part_name}  ({len(loaded)}/{len(modules)})")

        # ─── Загруженные модули ───
        for mod in loaded:
            print(f"      ✅ {mod}.py")

        # ─── Не загруженные модули ───
        for mod in not_loaded:
            print(f"      ⚠️  {mod}.py  (не загружен)")

        print()

    # ─── Итог ───
    print("─" * 80)
    print(f"  Загружено: {total_loaded} / {total_modules}")
    print("─" * 80)

# ============================================================================
# ЧАСТЬ 9c: ПОЛНЫЕ ОТЧЁТЫ МОДУЛЕЙ
# ============================================================================
#
# ЗАЧЕМ:
#   Вызвать print_report() (или аналогичную функцию) у каждого модуля,
#   чтобы в одном прогоне main_unified.py увидеть ПОЛНЫЙ вывод
#   каждого модуля — как будто запустили `python module.py`.
#
# КАК РАБОТАЕТ:
#   Для каждого модуля указана:
#     • функция отчёта (например, print_report, run_full_report, main)
#     • имя для отображения
#   Функция пытается импортировать модуль и вызвать указанную функцию.
#   Если модуль не загружен или функции нет — печатает предупреждение.
#
# ФОРМАТ ЗАПИСИ:
#   (module_name, func_name, title, kwargs, is_instance_method)
#
#   • module_name        — имя .py-файла без расширения
#   • func_name          — имя функции отчёта (может быть 'Class.method')
#   • title              — заголовок для печати
#   • kwargs             — словарь параметров функции отчёта
#   • is_instance_method — True, если это метод экземпляра класса
#                          (тогда создаётся экземпляр и вызывается метод)
#                          False — если это функция модуля или статический метод
#
# ВАЖНО:
#   • Некоторые модули НЕ имеют отдельной функции отчёта — пропускаем.
#   • Некоторые модули могут требовать параметры — передаём kwargs.
#   • Все ошибки ловятся и печатаются как ⚠️, не ломая прогон.
#
# ИСТОРИЯ:
#   • v1.3 — исходная версия (13 модулей).
#   • v1.4 — расширено до 38 модулей.                             ← ЭТА ВЕРСИЯ
#            Добавлена поддержка методов экземпляра (is_instance_method).
#            Добавлены: core, atlas, physics, statistical, cross_validation,
#            predictions, sensitivity, sensitivity2, benchmark,
#            topological_origin, extra_constants, system_metrics,
#            quantum_gravity, compactification, все research_*,
#            membrane_core, membrane_wave_function, genome_writer.
#            Исправлена ошибка elliptic_curve (было TypeError: missing 'self').
# ============================================================================

MODULE_REPORTS = [
    # ═══════════════════════════════════════════════════════════════════════
    # 1. ЯДРО (фундамент модели)
    # ═══════════════════════════════════════════════════════════════════════
    # core.py — константы, мастер-полином, Ходж, заряды CY₃
    ('core',                     'print_core_report',            'ЯДРО (core.py)',                {}, False),
    # atlas.py — атлас алгебр Ли, корневые системы, Φ⁺, связи E₆/E₇/E₈
    ('atlas',                    'print_atlas_report',           'АТЛАС АЛГЕБР ЛИ (atlas.py)',    {}, False),
    # physics.py — Юкава, GUT, Хиггс, космология, тетраэдр, петли, g-2
    ('physics',                  'print_physics_report',         'ФИЗИКА (physics.py)',           {}, False),

    # ═══════════════════════════════════════════════════════════════════════
    # 2. БАЗОВАЯ МОДЕЛЬ (main.py)
    # ═══════════════════════════════════════════════════════════════════════
    # Анализ вакуума SU(26) — три точки, вырождение, стабильность
    ('vacuum_analysis',          'print_report',                 'АНАЛИЗ ВАКУУМА',                {}, False),
    # Скан атласа — n_unique_cartan, классификация размерностей QNM
    ('atlas_scan',               'AtlasScan.print_report',       'СКАН АТЛАСА',                   {}, True),
    # Магические числа ядерной физики — 20, 28, 50, 82, 126
    ('magic_numbers',            'MagicNumbers.print_report',    'МАГИЧЕСКИЕ ЧИСЛА',              {}, True),
    # Структура кварков — n_q через n(k), связь b ↔ A(Pb)
    ('quark_structure',          'QuarkStructure.print_report',  'СТРУКТУРА КВАРКОВ',             {}, True),
    # Модулярные формы — ряды Эйзенштейна, коэффициенты Рамануджана
    ('modular_forms',            'ModularForms.print_report',    'МОДУЛЯРНЫЕ ФОРМЫ',              {}, True),
    # Эллиптическая кривая и BSD — L', bad prime 617, тёмная материя
    # ВАЖНО: is_instance_method=True, т.к. print_report — метод экземпляра
    ('elliptic_curve',           'EllipticCurveBSD.print_report','ЭЛЛИПТИЧЕСКАЯ КРИВАЯ',          {}, True),
    # Статистический анализ — Монте-Карло, Фишер, σ
    ('statistical',              'print_report',                 'СТАТИСТИКА',                    {}, False),
    # Кросс-валидация — LOO, K-Fold, Holdout
    ('cross_validation',         'print_report',                 'КРОСС-ВАЛИДАЦИЯ',               {}, False),
    # Предсказания — 30 предсказаний, 18 подтверждено
    ('predictions',              'print_report',                 'ПРЕДСКАЗАНИЯ',                  {}, False),
    # Чувствительность v2 
    ('sensitivity',              'print_report',                 'ЧУВСТВИТЕЛЬНОСТЬ v1',           {}, False),
    # Benchmark — сравнение с другими моделями
    ('benchmark',                'print_report',                 'BENCHMARK',                     {}, False),
    # Топологическое происхождение констант — CY₃ → β₁₂ → c₂ → B_i
    ('topological_origin',       'print_report',                 'ТОПОЛОГИЯ КОНСТАНТ',            {}, False),
    # Дополнительные константы
    ('extra_constants',          'print_report',                 'ДОП. КОНСТАНТЫ',                {}, False),
    # Системные метрики
    ('system_metrics',           'print_report',                 'СИСТЕМНЫЕ МЕТРИКИ',             {}, False),
    # Квантовая гравитация — M_Pl, G, Λ из SU(26)
    ('quantum_gravity',          'print_report',                 'КВАНТОВАЯ ГРАВИТАЦИЯ',          {}, False),
    # Компактификация 26D → 4D — разложение 675, массовая матрица
    ('compactification',         'print_report',                 'КОМПАКТИФИКАЦИЯ 26D → 4D',      {}, False),

    # ═══════════════════════════════════════════════════════════════════════
    # 3. ДОП-ПРОЕКТ (main_research.py)
    # ═══════════════════════════════════════════════════════════════════════
    # Топология, порог N_eff=8, формулы для 7
    ('research_topology',        'print_report',                 'RESEARCH: ТОПОЛОГИЯ',            {}, False),
    # Мембрана — гипотеза, спектр балка
    ('research_membrane',        'print_report',                 'RESEARCH: МЕМБРАНА',            {}, False),
    # Балк SU(58) — N_bulk=57, N_hidden=50, вложение
    ('research_bulk',            'print_report',                 'RESEARCH: БАЛК',                {}, False),
    # M-теория — 11D, E₈×E₈, дуальности
    ('research_m_theory',        'print_report',                 'RESEARCH: M-ТЕОРИЯ',            {}, False),
    # Ω⁴ — завершение, поля, действие, уравнения
    ('research_omega4',          'print_report',                 'RESEARCH: Ω⁴',                  {}, False),
    # Ω⁵⁸ — определение, спектр, гравитация, компактификация
    ('research_omega58',         'print_report',                 'RESEARCH: Ω⁵⁸',                 {}, False),
    # Ω-иерархия — 8/8 уровней
    ('research_hierarchy',       'print_report',                 'RESEARCH: ИЕРАРХИЯ',            {}, False),
    # Предсказания — экспериментальные, тёмная материя
    ('research_predictions',     'print_report',                 'RESEARCH: ПРЕДСКАЗАНИЯ',        {}, False),
    # Финальный отчёт доп-проекта
    ('research_final',           'print_report',                 'RESEARCH: ФИНАЛ',               {}, False),

    # ═══════════════════════════════════════════════════════════════════════
    # 4. FRG ↔ RGE
    # ═══════════════════════════════════════════════════════════════════════
    # Граница FRG — N_eff = 7 (между 7 и 8), инвариант g*+θ₁=4
    ('check_frg_limit',          'main',                         'FRG-ГРАНИЦА',                   {}, False),
    # FRG ↔ RGE — spread = 20.013016, STRUCTURAL_PROVEN
    ('omega4_frg_rge_unified',   'print_report',                 'FRG ↔ RGE',                     {}, False),

    # ═══════════════════════════════════════════════════════════════════════
    # 5. SUSY-СЕКТОР
    # ═══════════════════════════════════════════════════════════════════════
    # Факторизация Дарбу — H = A†A, изоспектральность
    ('darboux_susy',             'print_report',                 'ДАРБУ SUSY',                    {}, False),
    # SUSY Ω⁴ — M_SUSY = (2a−1)·M_EW = 3.69 ТэВ
    ('susy_omega4',              'print_report',                 'SUSY Ω⁴',                       {}, False),
    # Уникальность spread — D = 2^N_eff = 256
    ('spread_uniqueness',        'print_report',                 'УНИКАЛЬНОСТЬ SPREAD',           {}, False),
    # Вывод M_SUSY — спектр Пёшля-Теллера, ΔE = 2a−1
    ('msusy_derivation',         'run_full_report',              'ВЫВОД M_SUSY',                  {}, False),

    # ═══════════════════════════════════════════════════════════════════════
    # 6. МЕМБРАНА
    # ═══════════════════════════════════════════════════════════════════════
    # Ядро мембраны — 1D-модель, SUSY QM, 4D-базис, 28 констант
    ('membrane_core',            'print_report',                 'МЕМБРАНА: ЯДРО',                {}, False),
    # Волновая функция мембраны — CLI, полный прогон
    ('membrane_wave_function',   'print_report',                 'МЕМБРАНА: ВОЛНОВАЯ ФУНКЦИЯ',    {}, False),
    # Связь вакуума и мембраны — три точки ↔ три режима a
    ('membrane_vacuum_link',     'print_vacuum_membrane_report', 'СВЯЗЬ ВАКУУМА И МЕМБРАНЫ',      {}, False),

    # ═══════════════════════════════════════════════════════════════════════
    # 7. ГЕНОМ
    # ═══════════════════════════════════════════════════════════════════════
    # Генератор автогенома — список модулей для загрузки в REGISTRY
    ('genome_writer',            'print_modules_list',           'ГЕНОМ: СПИСОК МОДУЛЕЙ',         {}, False),
]

def print_all_module_reports(verbose: bool = True) -> dict:
    """
    Вызвать print_report() у каждого модуля из MODULE_REPORTS.

    ПОДДЕРЖКА МЕТОДОВ ЭКЗЕМПЛЯРА:
        Если is_instance_method=True, то:
          1. Имя модуля импортируется.
          2. Класс берётся из func_name (например 'EllipticCurveBSD.print_report').
          3. Создаётся экземпляр класса.
          4. Вызывается метод экземпляра.

        Если is_instance_method=False, но func_name содержит точку
        ('Class.method'), то метод вызывается у САМОГО КЛАССА
        (это работает для @staticmethod и @classmethod).

    ПОДДЕРЖКА СТАРОГО ФОРМАТА:
        Если запись содержит 4 элемента (без is_instance_method),
        считается, что is_instance_method=False.

    ВОЗВРАЩАЕТ:
        dict со статистикой:
          • 'called'      — список модулей, чьи функции вызваны успешно
          • 'failed'      — список модулей с ошибками (ImportError, Exception)
          • 'skipped'     — список модулей без нужной функции (AttributeError)
          • 'n_called'    — число успешных
          • 'n_total'     — общее число
    """
    if not verbose:
        return {'called': [], 'failed': [], 'skipped': [],
                'n_called': 0, 'n_total': 0}

    print()
    print("=" * 80)
    print("  📋 ПОЛНЫЕ ОТЧЁТЫ МОДУЛЕЙ")
    print("=" * 80)
    print()
    print(f"  Всего модулей с отчётами: {len(MODULE_REPORTS)}")
    print()

    called = []
    failed = []
    skipped = []

    for i, entry in enumerate(MODULE_REPORTS, 1):
        # ─── Поддержка старого (4) и нового (5) формата записи ───
        if len(entry) == 4:
            mod_name, func_name, title, kwargs = entry
            is_instance = False
        else:
            mod_name, func_name, title, kwargs, is_instance = entry

        # ─── Заголовок отчёта ───
        print()
        print("─" * 80)
        print(f"  [{i}/{len(MODULE_REPORTS)}] ОТЧЁТ: {title}")
        print(f"          Модуль:  {mod_name}.py")
        print(f"          Функция: {func_name}()")
        if is_instance:
            print(f"          Тип:     метод экземпляра")
        print("─" * 80)
        print()

        try:
            # ─── Импорт модуля ───
            import importlib
            mod = importlib.import_module(mod_name)

            # ─── Разбор func_name ───
            parts = func_name.split('.')

            if len(parts) == 1:
                # Простая функция модуля: module.print_report()
                func = getattr(mod, parts[0])
                func(**kwargs)

            elif len(parts) == 2:
                # Класс.метод: module.Class.method()
                cls = getattr(mod, parts[0])

                if is_instance:
                    # Создаём экземпляр и вызываем метод у него
                    # (для обычных методов, требующих self)
                    instance = cls()
                    method = getattr(instance, parts[1])
                else:
                    # Вызываем метод у класса
                    # (для @staticmethod и @classmethod)
                    method = getattr(cls, parts[1])

                method(**kwargs)

            else:
                raise AttributeError(
                    f"Неподдерживаемый формат func_name: {func_name}"
                )

            called.append(mod_name)
            print(f"\n  ✅ {mod_name}.{func_name}() — OK")

        except ImportError as e:
            # Модуль не найден или сломан импорт
            failed.append((mod_name, f'ImportError: {e}'))
            print(f"\n  ❌ {mod_name}.py — ImportError: {e}")

        except AttributeError as e:
            # Модуль загружен, но нужной функции/класса/метода нет
            # Это не критично — просто пропускаем
            skipped.append((mod_name, str(e)))
            print(f"\n  ⚠️  {mod_name}.{func_name}() — пропущен: {e}")

        except Exception as e:
            # Любая другая ошибка при вызове функции отчёта
            failed.append((mod_name, f'{type(e).__name__}: {e}'))
            print(f"\n  ❌ {mod_name}.{func_name}() — "
                  f'{type(e).__name__}: {e}')

    # ─── Итог ───
    print()
    print("=" * 80)
    print("  📋 ИТОГ ПО ОТЧЁТАМ МОДУЛЕЙ")
    print("=" * 80)
    print()
    print(f"  Успешно:   {len(called)} / {len(MODULE_REPORTS)}")

    if called:
        print(f"  Вызваны:")
        for mod in called:
            print(f"    ✅ {mod}.py")

    if skipped:
        print(f"\n  Пропущены (нет функции): {len(skipped)}")
        for mod, err in skipped:
            print(f"    ⚠️  {mod}.py — {err}")

    if failed:
        print(f"\n  Ошибки:    {len(failed)}")
        for mod, err in failed:
            print(f"    ❌ {mod}.py — {err}")

    print("=" * 80)

    return {
        'called': called,
        'failed': failed,
        'skipped': skipped,
        'n_called': len(called),
        'n_total': len(MODULE_REPORTS),
    }

# ============================================================================
# ЧАСТЬ 13: АУДИТ КОНФЛИКТОВ МЕЖДУ МОДУЛЯМИ
# ============================================================================
#
# ЗАЧЕМ:
#   Найти константы, которые в разных модулях имеют разные значения.
#   Авторитетный источник — core.REGISTRY.
#
# ============================================================================

def audit_conflicts(
    unified: Dict[str, Any],
    threshold: float = 0.1,
    max_show: int = 100,
) -> Dict[str, Any]:
    """
    Поиск конфликтов между модулями: где значения расходятся > threshold%.
    """
    print()
    print("=" * 80)
    print("  🔍 АУДИТ КОНФЛИКТОВ МЕЖДУ МОДУЛЯМИ")
    print("=" * 80)
    print()
    print(f"  Порог расхождения: {threshold}%")
    print()

    # ─── 1. Загружаем REGISTRY ─────────────────────────────────────────
    try:
        from core import REGISTRY
    except ImportError as e:
        print(f"  ⚠️  REGISTRY недоступен: {e}")
        print("=" * 80)
        return {'available': False, 'error': str(e)}

    registry = REGISTRY.all()

    if not registry:
        print("  ⚠️  REGISTRY пустой")
        print("=" * 80)
        return {'available': False, 'error': 'empty_registry'}

    print(f"  Констант в REGISTRY: {len(registry)}")
    print()

    # ─── 2. Собираем числовые значения из unified ──────────────────────
    all_values: Dict[str, list] = {}
    n_scanned = [0]

    # ─── КЛЮЧИ-ИСКЛЮЧЕНИЯ: семейство Ω_DM ───
    # Три ключа — три разных определения Ω_DM, НЕ конфликты:
    #   ec_omega_DM_TOPOLOGICAL — топологическое (elliptic_curve.py)
    #   DMP_OMEGA_DM_MEMBRANE  — мембранное (research_predictions.py)
    #   OMEGA_DM_CONSISTENT    — согласованное (среднее двух)
    #
    # ВАЖНО: key может быть int (например, {0: 1, 13: 4} для 4D-уровней).
    # Проверяем только строковые ключи.
    OMEGA_DM_KEYS = (
        'ec_omega_DM_TOPOLOGICAL',
        'DMP_OMEGA_DM_MEMBRANE',
        'OMEGA_DM_CONSISTENT',
    )

    def _is_excluded(key) -> bool:
        """True, если ключ относится к семейству Ω_DM (не конфликт)."""
        return isinstance(key, str) and any(
            k in key for k in OMEGA_DM_KEYS
        )

    def scan(d, path: str = '') -> None:
        if not isinstance(d, dict):
            return
        for key, val in d.items():
            # Пропускаем служебные ключи (начинаются с '_')
            if isinstance(key, str) and key.startswith('_'):
                continue

            # Пропускаем семейство Ω_DM (три определения — не конфликты)
            if _is_excluded(key):
                continue

            full_path = f'{path}.{key}' if path else str(key)

            if isinstance(val, bool):
                continue

            if isinstance(val, (int, float)):
                n_scanned[0] += 1
                all_values.setdefault(str(key), []).append((val, full_path))
            elif isinstance(val, dict):
                scan(val, full_path)
            elif isinstance(val, list):
                n_scanned[0] += 1
                all_values.setdefault(str(key) + '.__len__', []).append(
                    (len(val), full_path)
                )

    scan(unified, 'unified')

    print(f"  Числовых значений в unified: {n_scanned[0]}")
    print(f"  Уникальных имён ключей:      {len(all_values)}")
    print()
    # ─── 3. Сравниваем с REGISTRY ──────────────────────────────────────
    conflicts = []
    n_compared = 0

    for name, data in registry.items():
        if name not in all_values:
            continue
        ref_value = data.get('value')
        if not isinstance(ref_value, (int, float)):
            continue
        if isinstance(ref_value, bool):
            continue
        if ref_value == 0:
            continue

        for value, source in all_values[name]:
            if not isinstance(value, (int, float)):
                continue
            if isinstance(value, bool):
                continue
            n_compared += 1
            rel_diff = abs(value - ref_value) / abs(ref_value) * 100
            if rel_diff > threshold:
                conflicts.append({
                    'name': name,
                    'registry_value': ref_value,
                    'other_value': value,
                    'source': source,
                    'rel_diff_percent': rel_diff,
                    'registry_source': data.get('source', '?'),
                })

    # ─── 4. Печатаем ───────────────────────────────────────────────────
    print("─" * 80)

    if not conflicts:
        print("  ✅ КОНФЛИКТОВ НЕ НАЙДЕНО")
        print()
        print(f"  Проверено значений: {n_compared}")
        print("=" * 80)
        return {
            'available': True,
            'n_conflicts': 0,
            'conflicts': [],
            'n_checked': n_compared,
        }

    conflicts.sort(key=lambda x: -x['rel_diff_percent'])

    print(f"  ⚠️  НАЙДЕНО КОНФЛИКТОВ: {len(conflicts)}")
    print()

    print(f"  {'#':>3} | {'Константа':<22} | {'REGISTRY':>15} | "
          f"{'Другой':>15} | {'Откл.%':>9}")
    print(f"  {'-'*3} | {'-'*22} | {'-'*15} | {'-'*15} | {'-'*9}")

    for i, c in enumerate(conflicts[:max_show], 1):
        name_disp = c['name'][:22]

        def fmt_val(v):
            if v == 0:
                return '0'
            av = abs(v)
            if av >= 1e6 or av < 1e-3:
                return f'{v:.3e}'
            return f'{v:.6g}'

        print(f"  {i:>3} | {name_disp:<22} | "
              f"{fmt_val(c['registry_value']):>15} | "
              f"{fmt_val(c['other_value']):>15} | "
              f"{c['rel_diff_percent']:>8.4f}%")

    if len(conflicts) > max_show:
        print(f"  ... и ещё {len(conflicts) - max_show}")

    print()

    print("─" * 80)
    print("  ── ИСТОЧНИКИ КОНФЛИКТОВ (топ-20) ──")
    print()

    for i, c in enumerate(conflicts[:20], 1):
        print(f"  [{i}] {c['name']}")
        print(f"      REGISTRY: {c['registry_value']:.6g}  "
              f"(источник: {c['registry_source']})")
        print(f"      Другой:   {c['other_value']:.6g}  "
              f"(путь: {c['source']})")
        print(f"      Отклонение: {c['rel_diff_percent']:.4f}%")
        print()

    if len(conflicts) > 20:
        print(f"  ... и ещё {len(conflicts) - 20} конфликтов")
        print()

    print("─" * 80)
    print("  ── ИТОГ ──")
    print()
    print(f"  Констант в REGISTRY:        {len(registry)}")
    print(f"  Проверено значений:         {n_compared}")
    print(f"  Конфликтов найдено:         {len(conflicts)}")

    n_critical = sum(1 for c in conflicts if c['rel_diff_percent'] > 10)
    n_major = sum(1 for c in conflicts if 1 <= c['rel_diff_percent'] <= 10)
    n_minor = sum(1 for c in conflicts if c['rel_diff_percent'] < 1)

    print(f"    • Критических (> 10%):    {n_critical}")
    print(f"    • Значимых (1–10%):       {n_major}")
    print(f"    • Мелких (< 1%):          {n_minor}")
    print()

    print("  РЕКОМЕНДАЦИЯ:")
    if n_critical > 0:
        print("    🔴 СНАЧАЛА согласовать критические конфликты.")
        print("       Авторитетный источник — core.REGISTRY.")
    elif n_major > 0:
        print("    🟡 Согласовать значимые конфликты.")
    elif n_minor > 0:
        print("    🟢 Только мелкие конфликты — можно игнорировать.")
    else:
        print("    ✅ Всё согласовано.")

    print()
    print("=" * 80)

    return {
        'available': True,
        'n_conflicts': len(conflicts),
        'conflicts': conflicts,
        'n_checked': n_compared,
        'n_critical': n_critical,
        'n_major': n_major,
        'n_minor': n_minor,
    }

def main():
    """Полный объединённый прогон."""
    args = sys.argv[1:]

    save_json = '--no-save' not in args
    no_genome = '--no-genome' in args
    verbose = '--quiet' not in args

    print("=" * 80)
    print("  Ω⁴-UNIFIED — ОБЪЕДИНЁННЫЙ MAIN (main_unified.py) — v1.3")
    print("=" * 80)
    print()
    print("  Шесть частей:")
    print("    1. Базовая модель       (main.py)                     — 164/164")
    print("    2. Доп-проект           (main_research.py v3.0)       — 9 файлов")
    print("    3. FRG-диагностика      (check_frg_limit.py v2.1)     — граница 7/8")
    print("    4. FRG ↔ RGE анализ     (omega4_frg_rge_unified.py)   — spread = 20.013")
    print("    5. SUSY-сектор          (darboux + susy_omega4 + spread_uniqueness)")
    print("    6. Волновая функция     (membrane_core + membrane_wave_function")
    print("                             + membrane_vacuum_link)")
    print()
    print("  Порядок:")
    print("    1 → 2 → 3 → 4 → 5 → 6 → unified_results.json → genome_writer")
    print()
    print("  Аргументы:")
    print("    --no-genome  — не генерировать автогеном")
    print("    --no-save    — не сохранять unified_results.json")
    print("    --quiet      — минимум вывода")
    print()
    print("=" * 80)

    # ─── 1. Базовая модель ─────────────────────────────────────────────
    base_results = run_base_model(verbose=verbose)

    # ← печатаем все 174 проверки
    if verbose:
        print_base_checks_table(base_results)

    # ─── 2. Доп-проект ─────────────────────────────────────────────────
    research_data = run_research_project(verbose=verbose)

    # ─── 3. FRG-диагностика ────────────────────────────────────────────
    frg_boundary = run_frg_diagnostics(verbose=verbose)

    # ─── 4. FRG ↔ RGE АНАЛИЗ ───────────────────────────────────────────
    frg_rge_data = run_frg_rge_analysis(verbose=verbose)

    # ─── 5. SUSY-СЕКТОР ────────────────────────────────────────────────
    susy_data = run_susy_sector(verbose=verbose)

    # ─── 6. ★ ВОЛНОВАЯ ФУНКЦИЯ МЕМБРАНЫ ────────────────────────────────
    membrane_data = run_membrane_sector(verbose=verbose)

    # ─── 7. Объединённый JSON ──────────────────────────────────────────
    #
    # Собираем unified_results.json из результатов всех шести частей.
    #
    # ВАЖНО:
    #   • unified = None — инициализация ДО if save_json:
    #     (нужна, чтобы переменная unified существовала,
    #      даже если save_json=False — для аудита конфликтов в конце)
    #   • build_unified_json() — собирает словарь
    #   • __summary__ (служебный ключ) — удаляем из копии,
    #     чтобы он не попал в финальный JSON
    #
    unified = None
    if save_json:
        unified = build_unified_json(
            base_results=base_results,
            research_data=research_data,
            frg_boundary=frg_boundary,
            frg_rge_data=frg_rge_data,
            susy_data=susy_data,
            membrane_data=membrane_data,
        )
        # Убираем служебный ключ __summary__ из копии base_model.results
        results_copy = dict(unified['base_model']['results'])
        results_copy.pop('__summary__', None)
        unified['base_model']['results'] = results_copy

        # Сохраняем unified_results.json на диск
        save_unified_json(
            unified,
            output_path='unified_results.json',
            verbose=verbose,
        )

    # ─── 8. Финальная сводка ───────────────────────────────────────────
    print_final_summary(
        base_results, research_data, frg_boundary, frg_rge_data,
        susy_data=susy_data,
        membrane_data=membrane_data,
    )

    # ─── 9. Автогеном (один раз, в самом конце) ────────────────────────
    if not no_genome:
        try:
            from genome_writer import write_genome
            from core import REGISTRY

            print()
            print("─" * 80)
            print("  🌱 ГЕНЕРАЦИЯ АВТОМАТИЧЕСКОГО ГЕНОМА")
            print("─" * 80)
            print(f"  Источник: core.REGISTRY  ({len(REGISTRY)} констант)")
            write_genome(
                output_path='constants_genome_omega4.py',
                dry_run=False,
                verbose=True,
            )
        except ImportError as e:
            print(f"\n  ⚠️ genome_writer.py недоступен: {e}")
        except Exception as e:
            print(f"\n  ⚠️ Ошибка генерации автогенома: {e}")

    # ─── 10. Итоговый список загруженных модулей ───────────────────────
    #
    # Показывает, какие Python-модули проекта реально загрузились.
    # Вызывается ТОЛЬКО если не задан флаг --quiet.
    #
    if verbose:
        print_loaded_modules_summary()

    # ─── 11. Полные отчёты модулей ─────────────────────────────────────
    #
    # Вызывает print_report() (или аналогичную функцию) каждого модуля,
    # чтобы в одном прогоне увидеть ПОЛНЫЙ вывод каждого модуля.
    #
    # ВАЖНО: показывается только при --verbose (не при --quiet).
    #
    if verbose:
        print_all_module_reports(verbose=True)

    # ─── 12. Аудит конфликтов между модулями ───────────────────────────
    #
    # Сравнивает значения из core.REGISTRY (авторитет) со значениями,
    # собранными в unified_results.json (данные из 39 модулей).
    #
    # Находит константы, которые в разных модулях имеют разные значения
    # (например, M_26 = 3467 в core.py, но 5270 в physics.py).
    #
    # ВЫЗЫВАЕТСЯ:
    #   • только при verbose (не при --quiet)
    #   • только если unified_results.json был собран (save_json=True)
    #
    # ВОЗВРАЩАЕТ:
    #   dict с числом конфликтов и их деталями (для дальнейшего анализа).
    #
    if verbose and unified is not None:
        audit_conflicts(unified, threshold=0.1)

    # ─── 13. Итоговое сообщение ────────────────────────────────────────
    #
    # Сообщает, что прогон полностью завершён.
    #
    if verbose:
        print()
        print("=" * 80)
        print("  ✅ main_unified.py — ПРОГОН ПОЛНОСТЬЮ ЗАВЕРШЁН")
        print("=" * 80)
        print()


if __name__ == '__main__':
    main()
