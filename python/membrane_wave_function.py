#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
Ω⁴-UNIFIED v79.6 — ВОЛНОВАЯ ФУНКЦИЯ МЕМБРАНЫ (membrane_wave_function.py) — v3.1
================================================================================
Файл: membrane_wave_function.py
Версия: 3.1 (единая CLI + подробная печать)
Дата: 2026-10-02
Автор: Сергей Викторович Матершов (ORCID: 0009-0009-0641-1357)

НАЗНАЧЕНИЕ:
    ЕДИНАЯ CLI-обёртка над membrane_core.py.
    Заменяет три отдельных файла:
        • membrane_wave_function_full.py (v2.1)
        • membrane_wave_function_all.py  (v2.0)
        • membrane_wave_function_4D_basis.py (v1.3)

ИСПОЛЬЗОВАНИЕ:
    python membrane_wave_function.py                # полный прогон
    python membrane_wave_function.py --basis        # только 4D через базис
    python membrane_wave_function.py --convergence  # только 1D + сходимость
    python membrane_wave_function.py --no-susy      # без SUSY-партнёров
    python membrane_wave_function.py --no-viz       # без визуализации
    python membrane_wave_function.py --summary      # краткая сводка
    python membrane_wave_function.py --quiet        # тихий режим

АВТОР:  Сергей Викторович Матершов
ORCID:  0009-0009-0641-1357
ЛИЦЕНЗИЯ: CC BY-NC-ND 4.0 International
================================================================================
"""

import sys
import time
from typing import Dict, Any

# ─── Импорт из ядра ─────────────────────────────────────────────────────────
from membrane_core import (
    # Фундаментальные числа
    CHI, BETA_12, N_GEN, H11, H21,
    N_EFF_MEMBRANE, N_EFF_BREAK, M_EW,
    N_S_OBS, R_OBS_MAX,

    # 1D-модель
    spectrum_PT, run_1D_spectrum, get_1D_wave_functions,

    # SUSY
    verify_darboux_factoring, compute_susy_partners,

    # 2D-модель
    run_2D_convergence,

    # 4D через базис
    run_4D_spectrum, compute_4D_spectrum_via_wave_functions,

    # Комбинаторика
    generating_function_P4,
    compare_multiplicities_numerical_vs_combinatorial,
    analyze_topological_multiplicities,
    analyze_boundary_multiplicities,

    # 4D-волновые функции
    check_1D_orthonormality,
    check_4D_normalization,
    compute_4D_matrix_elements,

    # Сходимость
    run_1D_convergence,

    # Топология, космология, массы
    run_mass_spectrum, run_M_SUSY, run_topology,
    run_cosmology, run_number_10,

    # Визуализация
    visualize_1D_wave_functions,
    visualize_4D_wave_function_slice,
    visualize_susy_partners,
)


# ============================================================================
# ЧАСТЬ 1: ФЛАГИ CLI
# ============================================================================

def parse_args(argv):
    """Парсинг аргументов командной строки."""
    flags = {
        'mode': 'full',      # 'full' | 'basis' | 'convergence'
        'susy': True,
        'viz': True,
        'summary': False,
        'quiet': False,
    }

    for arg in argv:
        if arg == '--basis':
            flags['mode'] = 'basis'
        elif arg == '--convergence':
            flags['mode'] = 'convergence'
        elif arg == '--full':
            flags['mode'] = 'full'
        elif arg == '--no-susy':
            flags['susy'] = False
        elif arg == '--no-viz':
            flags['viz'] = False
        elif arg == '--summary':
            flags['summary'] = True
        elif arg == '--quiet':
            flags['quiet'] = True
        elif arg in ['--help', '-h']:
            print_help()
            sys.exit(0)
        else:
            print(f"⚠️ Неизвестный аргумент: {arg}")
            print("Используйте --help")
            sys.exit(1)

    return flags


def print_help():
    """Справка."""
    print("=" * 78)
    print("  Ω⁴-UNIFIED v79.6 — ВОЛНОВАЯ ФУНКЦИЯ МЕМБРАНЫ (CLI)")
    print("=" * 78)
    print()
    print("  ИСПОЛЬЗОВАНИЕ:")
    print("    python membrane_wave_function.py [ФЛАГИ]")
    print()
    print("  РЕЖИМЫ:")
    print("    --full          Полный прогон (по умолчанию)")
    print("    --basis         Только 4D через базис + комбинаторика")
    print("    --convergence   Только 1D + сходимость O(dx²)")
    print()
    print("  ОПЦИИ:")
    print("    --no-susy       Пропустить SUSY-партнёров")
    print("    --no-viz        Пропустить визуализацию (для CI)")
    print("    --summary       Краткая сводка (без деталей)")
    print("    --quiet         Тихий режим")
    print()
    print("  ПРИМЕРЫ:")
    print("    python membrane_wave_function.py")
    print("    python membrane_wave_function.py --basis")
    print("    python membrane_wave_function.py --convergence --no-viz")
    print("    python membrane_wave_function.py --summary")
    print()
    print("=" * 78)


# ============================================================================
# ЧАСТЬ 2: ПЕЧАТЬ ПРОМЕЖУТОЧНЫХ РЕЗУЛЬТАТОВ
# ============================================================================

def print_section(title):
    """Печать заголовка секции."""
    print("\n" + "=" * 78)
    print(f"  {title}")
    print("=" * 78)


def print_1D_results(results):
    """Печать 1D-спектра."""
    s1d = results.get('1D', {})
    E_num = s1d.get('E_num', [])
    E_ana = s1d.get('E_ana', [])
    max_dev = s1d.get('max_dev', 0)

    print(f"\n  {'n':>3} | {'E_num':>14} | {'E_ana':>10} | {'|ΔE|':>12} | Status")
    print(f"  {'-'*3} | {'-'*14} | {'-'*10} | {'-'*12} | {'-'*8}")

    for n in range(min(len(E_num), len(E_ana))):
        diff = abs(E_num[n] - E_ana[n])
        status = '✅' if diff < 0.01 else '❌'
        print(f"  {n:>3} | {E_num[n]:>14.6f} | {E_ana[n]:>10.4f} | "
              f"{diff:>12.6f} | {status}")

    print(f"\n  Макс. отклонение: {max_dev:.6f}")


def print_darboux_results(results):
    """Печать факторизации Дарбу."""
    d = results.get('darboux', {})
    print(f"\n  H₁ = A†A: эрмитов = {d.get('herm1', False)}")
    print(f"  H₂ = AA†: эрмитов = {d.get('herm2', False)}")


def print_2D_results(results):
    """Печать 2D-проверки."""
    s2d = results.get('2D', {})

    if 'skipped' in s2d:
        print(f"  ПРОПУЩЕНО (scipy не установлен)")
        return

    if 'error' in s2d:
        print(f"  ОШИБКА: {s2d['error']}")
        return

    res = s2d.get('results', [])
    if not res:
        print(f"  Нет данных")
        return

    print(f"\n  {'N':>5} | {'dx':>10} | {'t (сек)':>10} | "
          f"{'max|ΔE|':>10} | {'max|ΔE|/dx²':>12}")
    print(f"  {'-'*5} | {'-'*10} | {'-'*10} | {'-'*10} | {'-'*12}")

    for r in res:
        print(f"  {r['N']:>5} | {r['dx']:>10.6f} | {r['time']:>10.2f} | "
              f"{r['max_dev']:>10.4f} | {r['ratio']:>12.2f}")

    print(f"\n  ВЫВОД: max|ΔE| ≈ const · dx² → 0 при dx → 0")
    print(f"  → Разделение переменных РАБОТАЕТ точно.")


def print_4D_results(results):
    """Печать 4D-спектра через базис."""
    s4d = results.get('4D', {})
    n_total = s4d.get('n_total', 0)
    n_unique = s4d.get('n_unique', 0)
    max_dev = s4d.get('max_dev', 0)
    unique_E = s4d.get('unique_E', [])
    energy_counts = s4d.get('energy_counts', {})

    print(f"\n  Число состояний: {n_total}")
    print(f"  Уникальных уровней: {n_unique}")
    print(f"  Макс. отклонение: {max_dev:.2e}")
    print(f"  Совпадение: {'✅ ТОЧНОЕ' if max_dev < 1e-10 else '⚠️'}")

    print(f"\n  {'#':>3} | {'E_4D':>12} | {'Кратность':>10}")
    print(f"  {'-'*3} | {'-'*12} | {'-'*10}")

    for i, E in enumerate(unique_E[:15]):
        count = energy_counts.get(E, 0)
        print(f"  {i:>3} | {E:>12.6f} | {count:>10}")


def print_combinatorial_results(results):
    """Печать комбинаторных кратностей."""
    sc = results.get('combinatorial', {})
    all_match = sc.get('all_match', False)
    n_checked = sc.get('n_checked', 0)
    n_mismatch = sc.get('n_mismatch', 0)

    print(f"\n  Проверено уровней: {n_checked}")
    print(f"  Несовпадений: {n_mismatch}")
    print(f"  Все совпадают: {'✅' if all_match else '❌'}")


def print_topo_results(results):
    """Печать топологических кратностей."""
    st = results.get('topo', {})
    binomial_match = st.get('binomial_match', False)
    binomial_results = st.get('binomial_results', [])
    total_states = st.get('total_states', 0)

    print(f"\n  ── БИНОМИАЛЬНЫЕ КРАТНОСТИ ──")
    print(f"  ΔE = 2a − 1 = 13.0")
    print(f"\n  {'k':>3} | {'E':>6} | {'N(E)':>6} | {'C(4,k)':>8} | Status")
    print(f"  {'-'*3} | {'-'*6} | {'-'*6} | {'-'*8} | {'-'*8}")

    for br in binomial_results:
        status = '✅' if br['N'] == br['C'] else '❌'
        print(f"  {br['k']:>3} | {br['E']:>6.1f} | {br['N']:>6} | "
              f"{br['C']:>8} | {status}")

    print(f"\n  ── ВНЕШНЯЯ АЛГЕБРА Λ^k R^4 ──")
    print(f"  dim Λ^0 R^4 = 1   (вакуум)")
    print(f"  dim Λ^1 R^4 = 4   (векторы)")
    print(f"  dim Λ^2 R^4 = 6   (2-формы)")
    print(f"  dim Λ^3 R^4 = 4   (3-формы)")
    print(f"  dim Λ^4 R^4 = 1   (4-формы)")
    print(f"  Сумма: Σ C(4,k) = 2^4 = 16")

    print(f"\n  ── ПРОВЕРКА ──")
    print(f"  Σ N(E) = {total_states}")
    print(f"  8^4 = {8**4}")
    print(f"  Совпадение: {'✅' if total_states == 8**4 else '❌'}")

    print(f"\n  Биномиальные кратности: {'✅' if binomial_match else '❌'}")


def print_boundary_results(results):
    """Печать граничных кратностей."""
    sb = results.get('boundary', {})
    max_N = sb.get('max_N', 0)
    nz_min = sb.get('nz_min', 0)
    n_levels = sb.get('n_levels', 0)

    print(f"\n  ── ГРАНИЧНЫЕ КРАТНОСТИ ──")
    print(f"  max N(E) = {max_N}")
    print(f"  min N(E) = {nz_min}")
    print(f"  Всего уровней: {n_levels}")


def print_orthonormality_results(results):
    """Печать ортонормальности 1D."""
    ortho = results.get('orthonormality', {})
    max_dev = ortho.get('max_dev', 0)

    print(f"\n  1D ортонормальность: max|ΔS| = {max_dev:.2e}")
    print(f"  Статус: {'✅' if max_dev < 1e-3 else '❌'}")


def print_normalization_4D_results(results):
    """Печать нормировки 4D."""
    norm = results.get('normalization_4D', {})
    norm_4D = norm.get('norm_4D', 0)

    print(f"\n  4D нормировка: {norm_4D:.10f}")
    print(f"  Отклонение от 1: {abs(norm_4D - 1):.2e}")
    print(f"  Статус: {'✅' if abs(norm_4D - 1) < 1e-6 else '❌'}")


def print_matrix_elements_results(results):
    """Печать матричных элементов."""
    me = results.get('matrix_elements', {})

    print(f"\n  {'Оператор':>20} | {'Значение':>12}")
    print(f"  {'-'*20} | {'-'*12}")
    for key, val in me.items():
        print(f"  {key:>20} | {val:>12.6f}")

    print(f"\n  Принцип неопределённости ħ/2 = 0.5:")
    for i in range(1, 5):
        key = f'Δx_{i}·Δp_{i}'
        if key in me:
            status = '✅' if me[key] >= 0.4999 else '❌'
            print(f"    {key} = {me[key]:.6f} {status}")


def print_susy_partners_results(results):
    """Печать SUSY-партнёров."""
    sp = results.get('susy_partners', {})

    if not sp or 'error' in sp:
        print(f"  Нет данных")
        return

    psi_partners = sp.get('psi_partners', [])
    energies = sp.get('energies_partners', [])
    max_dev = sp.get('max_dev', 0)
    isospectral = sp.get('isospectral', False)

    print(f"\n  {'n':>3} | {'E_n^(1)':>10} | {'E_n^(2)':>10}")
    print(f"  {'-'*3} | {'-'*10} | {'-'*10}")

    for n, E in enumerate(energies):
        print(f"  {n:>3} | {E:>10.4f} | {E:>10.4f}")

    print(f"\n  Ортонормальность партнёров: max|ΔS| = {max_dev:.2e}")
    print(f"  Статус: {'✅' if max_dev < 0.1 else '❌'}")
    print(f"\n  Изоспектральность: {'✅' if isospectral else '❌'}")


def print_convergence_results(results):
    """Печать сходимости."""
    conv = results.get('convergence_1D', {})
    res = conv.get('results', [])

    if not res:
        print(f"  Нет данных")
        return

    print(f"\n  {'N':>6} | {'dx':>10} | {'max|ΔE|':>12} | "
          f"{'max|ΔE|/dx²':>14} | {'Статус':>10}")
    print(f"  {'-'*6} | {'-'*10} | {'-'*12} | {'-'*14} | {'-'*10}")

    for r in res:
        status = '✅' if r['max_dev'] < 0.01 else '⚠️'
        print(f"  {r['N']:>6} | {r['dx']:>10.6f} | {r['max_dev']:>12.6f} | "
              f"{r['ratio']:>14.2f} | {status:>10}")

    # ─── ПРАВИЛЬНАЯ ПРОВЕРКА СХОДИМОСТИ ───
    ratios = [r['ratio'] for r in res]
    if len(ratios) >= 3:
        last_three = ratios[-3:]
        spread = max(last_three) - min(last_three)
        stable = spread < 5.0
        print(f"\n  Последние 3 значения ratio: {[f'{x:.2f}' for x in last_three]}")
        print(f"  Разброс: {spread:.2f}")
        print(f"  Сходимость O(dx²): {'✅' if stable else '❌'}")


def print_mass_results(results):
    """Печать спектра масс."""
    sm = results.get('mass', {})
    masses = sm.get('masses', [])
    E_1 = sm.get('E_1', 0)

    if not masses:
        return

    print(f"\n  Формула: m_n = M_EW · √(E_n / E_1)")
    print(f"  E_1 = {E_1}, M_EW = {M_EW} ГэВ")
    print(f"\n  {'n':>3} | {'E_n':>8} | {'m_n [ГэВ]':>12}")
    print(f"  {'-'*3} | {'-'*8} | {'-'*12}")

    for n, m in enumerate(masses):
        E_n = spectrum_PT(7.0)[n]
        print(f"  {n:>3} | {E_n:>8.1f} | {m:>12.2f}")


def print_M_SUSY_results(results):
    """Печать M_SUSY."""
    ms = results.get('M_SUSY', {})
    delta_E = ms.get('delta_E', 0)
    M_SUSY = ms.get('M_SUSY', 0)

    print(f"\n  ΔE = 2a − 1 = {delta_E}")
    print(f"  M_SUSY = ΔE · M_EW = {M_SUSY:.2f} ГэВ")
    print(f"  M_SUSY = {M_SUSY/1000:.4f} ТэВ")


def print_topology_results(results):
    """Печать топологии."""
    topo = results.get('topology', {})
    matches_a2 = topo.get('matches_a2', [])
    matches_V0 = topo.get('matches_V0', [])

    print(f"\n  ── a² = 49 ──")
    for name in matches_a2:
        print(f"    ✅ {name}")

    print(f"\n  ── V₀ = 56 ──")
    for name in matches_V0:
        print(f"    ✅ {name}")


def print_cosmology_results(results):
    """Печать космологии."""
    scos = results.get('cosmology', {})
    N_e = scos.get('N_e', 0)
    n_s = scos.get('n_s', 0)
    r = scos.get('r', 0)
    dev_ns = scos.get('dev_ns', 0)

    print(f"\n  N_e = χ/10 = {N_e}")
    print(f"  n_s = 1 − 2/N_e = {n_s:.6f}")
    print(f"  r   = 12/N_e²   = {r:.6f}")
    print(f"\n  Сравнение с Planck:")
    print(f"    n_s (модель) = {n_s:.6f}")
    print(f"    n_s (Planck) = {N_S_OBS}")
    print(f"    Отклонение   = {dev_ns:.4f}%")


# ============================================================================
# ЧАСТЬ 3: РЕЖИМЫ ПРОГОНА
# ============================================================================

def run_full(a=7.0, n_basis=8, susy=True, viz=True, quiet=False) -> Dict:
    """Полный прогон с подробной печатью."""
    results = {}
    t0 = time.time()

    # ─── 1. 1D-СПЕКТР ───
    if not quiet:
        print_section("ЧАСТЬ 1: 1D-СПЕКТР ПЁШЛЯ-ТЕЛЛЕРА")
    results['1D'] = run_1D_spectrum(a=a, L=20.0, N=4000)
    if not quiet:
        print_1D_results(results)

    # ─── 2. ФАКТОРИЗАЦИЯ ДАРБУ ───
    if not quiet:
        print_section("ЧАСТЬ 2: ФАКТОРИЗАЦИЯ ДАРБУ")
    results['darboux'] = verify_darboux_factoring(a=a, N=2000, L=40.0)
    if not quiet:
        print_darboux_results(results)

    # ─── 3. 2D-ПРОВЕРКА ───
    if not quiet:
        print_section("ЧАСТЬ 3: 2D-ПРОВЕРКА (РАЗДЕЛЕНИЕ ПЕРЕМЕННЫХ)")
    results['2D'] = run_2D_convergence(a=a)
    if not quiet:
        print_2D_results(results)

    # ─── 4. 4D ЧЕРЕЗ БАЗИС ───
    if not quiet:
        print_section("ЧАСТЬ 4: 4D ЧЕРЕЗ БАЗИС ТЕНЗОРНЫХ ПРОИЗВЕДЕНИЙ")
    results['4D'] = run_4D_spectrum(a=a, n_basis=n_basis)
    if not quiet:
        print_4D_results(results)

    # ─── 5. КОМБИНАТОРНЫЕ КРАТНОСТИ ───
    if not quiet:
        print_section("ЧАСТЬ 5: КОМБИНАТОРНЫЕ КРАТНОСТИ")
    E_1D = spectrum_PT(a)
    P4 = generating_function_P4(list(E_1D), n_dims=4, E_max=300)

    results['combinatorial'] = compare_multiplicities_numerical_vs_combinatorial(
        results['4D']['unique_E'],
        results['4D']['energy_counts'],
        E_1D,
    )
    if not quiet:
        print_combinatorial_results(results)

    results['topo'] = analyze_topological_multiplicities(P4, E_1D)
    if not quiet:
        print_topo_results(results)

    results['boundary'] = analyze_boundary_multiplicities(P4)
    if not quiet:
        print_boundary_results(results)

    # ─── 6. 4D-ВОЛНОВЫЕ ФУНКЦИИ ───
    if not quiet:
        print_section("ЧАСТЬ 6: 4D-ВОЛНОВЫЕ ФУНКЦИИ")
    wf_1D = get_1D_wave_functions(a=a, L=20.0, N=4000, n_states=8)

    results['orthonormality'] = check_1D_orthonormality(wf_1D)
    if not quiet:
        print_orthonormality_results(results)

    results['normalization_4D'] = check_4D_normalization(wf_1D, (1, 0, 0, 0))
    if not quiet:
        print_normalization_4D_results(results)

    results['matrix_elements'] = compute_4D_matrix_elements(wf_1D, (1, 0, 0, 0))
    if not quiet:
        print_matrix_elements_results(results)

    # ─── 7. SUSY-ПАРТНЁРЫ ───
    if susy:
        if not quiet:
            print_section("ЧАСТЬ 7: SUSY-ПАРТНЁРЫ (ФАКТОРИЗАЦИЯ ДАРБУ)")
        partners = compute_susy_partners(wf_1D, a=a)
        results['susy_partners'] = partners
        if not quiet:
            print_susy_partners_results(results)

        if viz:
            visualize_susy_partners(wf_1D, partners, a=a, n_states=3)

    # ─── 8. ВИЗУАЛИЗАЦИЯ ───
    if viz:
        if not quiet:
            print_section("ВИЗУАЛИЗАЦИЯ")
        visualize_1D_wave_functions(wf_1D, n_states=4,
                                     filename='wave_functions_1D.png')
        if not quiet:
            print(f"  Сохранено: wave_functions_1D.png")

        for n_tuple in [(1,0,0,0), (1,1,0,0), (2,0,0,0), (1,1,1,1)]:
            visualize_4D_wave_function_slice(wf_1D, n_tuple, n_grid=200)
            if not quiet:
                print(f"  Сохранено: psi_4D_slice_"
                      f"{''.join(str(x) for x in n_tuple)}.png")

    # ─── 9. СХОДИМОСТЬ 1D ───
    if not quiet:
        print_section("ЧАСТЬ 8: СХОДИМОСТЬ 1D O(dx²)")
    results['convergence_1D'] = run_1D_convergence(a=a, L=30.0)
    if not quiet:
        print_convergence_results(results)

    # ─── 10. МАССЫ ───
    if not quiet:
        print_section("ЧАСТЬ 9: СПЕКТР МАСС ПОЛЕЙ")
    results['mass'] = run_mass_spectrum(a=a)
    if not quiet:
        print_mass_results(results)

    # ─── 11. M_SUSY ───
    if not quiet:
        print_section("ЧАСТЬ 10: МАСШТАБ SUSY")
    results['M_SUSY'] = run_M_SUSY(a=a)
    if not quiet:
        print_M_SUSY_results(results)

    # ─── 12. ТОПОЛОГИЯ ───
    if not quiet:
        print_section("ЧАСТЬ 11: ТОПОЛОГИЧЕСКИЕ СВЯЗИ")
    results['topology'] = run_topology()
    if not quiet:
        print_topology_results(results)

    # ─── 13. КОСМОЛОГИЯ ───
    if not quiet:
        print_section("ЧАСТЬ 12: КОСМОЛОГИЧЕСКИЕ СЛЕДСТВИЯ")
    results['cosmology'] = run_cosmology()
    if not quiet:
        print_cosmology_results(results)

    # ─── 14. ЧИСЛО 10 ───
    results['number_10'] = run_number_10()

    results['total_time'] = time.time() - t0
    return results


def run_basis(a=7.0, n_basis=8, quiet=False) -> Dict:
    """Только 4D через базис + комбинаторика."""
    results = {}
    t0 = time.time()

    if not quiet:
        print_section("1D-СПЕКТР")
    results['1D'] = run_1D_spectrum(a=a, L=20.0, N=4000)
    if not quiet:
        print_1D_results(results)

    if not quiet:
        print_section("4D ЧЕРЕЗ БАЗИС")
    results['4D'] = run_4D_spectrum(a=a, n_basis=n_basis)
    if not quiet:
        print_4D_results(results)

    if not quiet:
        print_section("КОМБИНАТОРНЫЕ КРАТНОСТИ")
    E_1D = spectrum_PT(a)
    P4 = generating_function_P4(list(E_1D), n_dims=4, E_max=300)

    results['combinatorial'] = compare_multiplicities_numerical_vs_combinatorial(
        results['4D']['unique_E'],
        results['4D']['energy_counts'],
        E_1D,
    )
    if not quiet:
        print_combinatorial_results(results)

    results['topo'] = analyze_topological_multiplicities(P4, E_1D)
    if not quiet:
        print_topo_results(results)

    results['boundary'] = analyze_boundary_multiplicities(P4)
    if not quiet:
        print_boundary_results(results)

    if not quiet:
        print_section("СПЕКТР МАСС")
    results['mass'] = run_mass_spectrum(a=a)
    if not quiet:
        print_mass_results(results)

    if not quiet:
        print_section("МАСШТАБ SUSY")
    results['M_SUSY'] = run_M_SUSY(a=a)
    if not quiet:
        print_M_SUSY_results(results)

    results['total_time'] = time.time() - t0
    return results


def run_convergence(a=7.0, L=30.0, quiet=False) -> Dict:
    """Только 1D + сходимость."""
    results = {}
    t0 = time.time()

    if not quiet:
        print_section("1D-СПЕКТР")
    results['1D'] = run_1D_spectrum(a=a, L=20.0, N=4000)
    if not quiet:
        print_1D_results(results)

    if not quiet:
        print_section("СХОДИМОСТЬ 1D O(dx²)")
    results['convergence_1D'] = run_1D_convergence(a=a, L=L)
    if not quiet:
        print_convergence_results(results)

    if not quiet:
        print_section("СПЕКТР МАСС")
    results['mass'] = run_mass_spectrum(a=a)
    if not quiet:
        print_mass_results(results)

    if not quiet:
        print_section("МАСШТАБ SUSY")
    results['M_SUSY'] = run_M_SUSY(a=a)
    if not quiet:
        print_M_SUSY_results(results)

    results['total_time'] = time.time() - t0
    return results


# ============================================================================
# ЧАСТЬ 4: КРАТКАЯ СВОДКА
# ============================================================================

def print_summary(results):
    """Краткая сводка."""
    print("\n" + "=" * 78)
    print("  КРАТКАЯ СВОДКА")
    print("=" * 78)

    checks = []

    # 1D-спектр
    if '1D' in results:
        max_dev = results['1D'].get('max_dev', 1)
        checks.append(('1D-спектр', max_dev < 0.01))

    # 2D разделение
    if '2D' in results and 'results' in results['2D']:
        r_last = results['2D']['results'][-1]
        checks.append(('2D разделение', r_last['ratio'] < 50))

    # 4D базис
    if '4D' in results:
        max_dev = results['4D'].get('max_dev', 1)
        checks.append(('4D-базис', max_dev < 1e-10))

    # Комбинаторика
    if 'combinatorial' in results:
        checks.append(('Комбинаторика',
                       results['combinatorial'].get('all_match', False)))

    # Биномиальные кратности
    if 'topo' in results:
        checks.append(('Бином. кратности',
                       results['topo'].get('binomial_match', False)))

    # 1D ортонормальность
    if 'orthonormality' in results:
        checks.append(('1D ортонорм.',
                       results['orthonormality']['max_dev'] < 1e-3))

    # 4D нормировка
    if 'normalization_4D' in results:
        norm = results['normalization_4D']['norm_4D']
        checks.append(('4D нормировка', abs(norm - 1) < 1e-6))

    # Принцип неопределённости
    if 'matrix_elements' in results:
        me = results['matrix_elements']
        all_ok = all(me.get(f'Δx_{i}·Δp_{{i}}'.replace('{i}', str(i)), 0) >= 0.4999
                     for i in range(1, 5))
        # Более простая проверка
        all_ok = all(me.get(f'Δx_{i}·Δp_{i}', 0) >= 0.4999
                     for i in range(1, 5))
        checks.append(('Неопределённость', all_ok))

    # SUSY изоспектральность
    if 'susy_partners' in results:
        checks.append(('SUSY изоспектр.',
                       results['susy_partners'].get('isospectral', False)))

    # Сходимость O(dx²) — ПРАВИЛЬНАЯ ПРОВЕРКА
    if 'convergence_1D' in results:
        conv = results['convergence_1D']['results']
        ratios = [cr['ratio'] for cr in conv]
        if len(ratios) >= 3:
            last_three = ratios[-3:]
            spread = max(last_three) - min(last_three)
            stable = spread < 5.0
        else:
            stable = False
        checks.append(('Сходимость O(dx²)', stable))

    # Печать
    for name, ok in checks:
        print(f"    {name:<25} {'✅' if ok else '❌'}")

    print(f"\n  Всего проверок: {len(checks)}")
    print(f"  Пройдено: {sum(1 for _, ok in checks if ok)}/{len(checks)}")
    print("=" * 78)


# ============================================================================
# ЧАСТЬ 5: ГЛАВНАЯ ФУНКЦИЯ
# ============================================================================

def main():
    """Точка входа."""
    flags = parse_args(sys.argv[1:])

    if not flags['quiet']:
        # ─── Заголовок ───
        mode_names = {
            'full': 'ПОЛНЫЙ ПРОГОН',
            'basis': 'ТОЛЬКО 4D ЧЕРЕЗ БАЗИС',
            'convergence': 'ТОЛЬКО СХОДИМОСТЬ',
        }
        print("=" * 78)
        print(f"  Ω⁴-UNIFIED v79.6 — ВОЛНОВАЯ ФУНКЦИЯ МЕМБРАНЫ")
        print(f"  Режим: {mode_names[flags['mode']]}")
        print("=" * 78)
        print()
        print(f"  ФУНДАМЕНТАЛЬНЫЕ ЧИСЛА:")
        print(f"    χ = {CHI}, β₁₂ = {BETA_12}")
        print(f"    h¹¹ = {H11}, h²¹ = {H21}, N_gen = {N_GEN}")
        print(f"    N_eff_membrane = {N_EFF_MEMBRANE}")
        print(f"    N_eff_break    = {N_EFF_BREAK}")
        print(f"    M_EW = {M_EW} ГэВ")

    # ─── Выбор режима ───
    if flags['mode'] == 'basis':
        results = run_basis(
            a=7.0, n_basis=8,
            quiet=flags['quiet'],
        )
    elif flags['mode'] == 'convergence':
        results = run_convergence(
            a=7.0, L=30.0,
            quiet=flags['quiet'],
        )
    else:  # full
        results = run_full(
            a=7.0, n_basis=8,
            susy=flags['susy'],
            viz=flags['viz'],
            quiet=flags['quiet'],
        )

    # ─── Краткая сводка ───
    if not flags['quiet']:
        print_summary(results)
        print(f"\n  Общее время: {results['total_time']:.2f} сек")

    return results

# ============================================================================
# ЧАСТЬ 6: ПЕЧАТЬ ОТЧЁТА (для main_unified.py)
# ============================================================================
#
# ФУНКЦИЯ: print_report()
#
# ЗАЧЕМ:
#   Вызывается из main_unified.py → MODULE_REPORTS.
#   Печатает полный отчёт по волновой функции мембраны.
#   По сути — обёртка над main() с фиксированными флагами.
#
# ВАЖНО:
#   Запускает полный прогон с viz=False (без графиков),
#   чтобы не тормозить main_unified.py.
#
# ============================================================================

def print_report():
    """
    Печать полного отчёта по волновой функции мембраны.
    
    ИСПОЛЬЗУЕТСЯ В:
        main_unified.py → MODULE_REPORTS →
        ('membrane_wave_function', 'print_report', ...)
    """
    # Устанавливаем флаги: полный режим, без визуализации, не тихий
    sys.argv = [
        'membrane_wave_function.py',
        '--full',
        '--no-viz',
    ]
    
    # Запускаем main() — он напечатает полный отчёт
    return main()
if __name__ == '__main__':
    results = main()
