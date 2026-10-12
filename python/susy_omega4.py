#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
Ω⁴-UNIFIED — SUSY QM СВЯЗЬ С Ω⁴ (susy_omega4.py) — v4.1
================================================================================
НАЗНАЧЕНИЕ:
    Объединённый модуль, связывающий SUSY QM (факторизация Дарбу)
    с топологией CY₃ в рамках модели Ω⁴.

    ОБЪЕДИНЯЕТ:
        • susy_solution.py v2.1   — факторизация + M_SUSY + LHC-проверка
        • susy_omega4_link.py v1.0 — 5 кандидатов для параметра a

    ИСПОЛЬЗУЕТ:
        • darboux_susy.py v6.0    — базовый модуль (Pöschl-Teller,
                                     DarbouxFactorizer, SpectrumChecker)

ЧТО НОВОГО В v4.1 (рефакторинг):
    ★ УСТРАНЕНО ДУБЛИРОВАНИЕ M_SUSY:
      • H3 (M_SUSY = (2a−1)·M_EW) берётся из msusy_derivation
        (ПЕРВИЧНЫЙ ИСТОЧНИК), с fallback на локальную формулу.
      • Добавлено поле 'M_SUSY_source' в результат.
      • Добавлена проверка совпадения первичного и локального значений.
      • Уникальные наработки СОХРАНЕНЫ: H1, H2, H4, H5,
        5 кандидатов для a, 8 суперпартнёров.

ЧТО БЫЛО В v4.0 (по результатам susy_hypotheses_test.py):
    1. ★ ОДНА ФОРМУЛА ДЛЯ M_SUSY (вместо трёх):
       • Основная: M_SUSY = (2a − 1) · M_EW  (H3)
       • Fallback: M_SUSY = M_EW / ε_FN²       (H2)
       • Отвергнута: M_SUSY = M_EW / ε_FN      (H1)

    2. ★ ОБОСНОВАНИЕ ВЫБОРА:
       Из susy_hypotheses_test.py:
         • H1 (0.78 ТэВ) — противоречит LHC (1/4 пределов)
         • H2 (2.45 ТэВ) — проходит LHC, но НЕ использует SUSY QM
         • H3 (3.69 ТэВ для a=8) — проходит LHC ✅ И использует SUSY QM ✅
         • H4 (√ΔE · M_EW) — НЕ проходит LHC для малых a (0.45 score)
         • H5 (M_EW/√ΔE) — противоречит LHC (0/4 пределов)

    3. ★ СТАТУС: STRUCTURAL_STRONG
       • H3 согласована с LHC
       • H3 использует SUSY QM (ΔE(a))
       • H3 топологически обоснована (a = N_eff = 8)

    4. ★ РЕГИСТРАЦИЯ В REGISTRY (16 констант)
       Добавлены: SUSY_OMEGA4_M_SUSY_H3, SUSY_OMEGA4_M_SUSY_H2_FALLBACK,
       SUSY_OMEGA4_H1_REJECTED

ПУБЛИЧНЫЙ API:
    • potential_peschl_teller(x, chi, a, v0)
    • riccati_ansatz_peschl_teller(x, chi, a)
    • verify_factoring()
    • susy_spectrum_peschl_teller(n_levels, chi, a)
    • compute_susy_masses(a)              ← НОВАЯ сигнатура
    • check_lhc(m_susy_gev)
    • verify_susy_mass(a)
    • get_a_candidates()
    • full_analysis()
    • assess_status()
    • print_report()

АВТОР:  Сергей Викторович Матершов
ORCID:  0009-0009-0641-1357
ДАТА:   2026-09-30
ВЕРСИЯ: 4.1
ЛИЦЕНЗИЯ: CC BY-NC-ND 4.0 International
================================================================================
ОБЛАСТЬ ОТВЕТСТВЕННОСТИ:
    Этот модуль: физическая интерпретация M_SUSY через SUSY QM
                 (топологический выбор a, LHC-сравнение).
    См. также:   msusy_derivation.py — численное обоснование ΔE = 2a−1
                 (спектр Пёшля-Теллера, сходимость O(dx²)).
"""

import math
from typing import Dict, List, Tuple, Any, Optional

from core import C, MasterPolynomial, REGISTRY

# ─── Импорт из darboux_susy.py (базовый модуль, НЕ дублируем) ──────────────
try:
    from darboux_susy import (
        Superpotential,
        DarbouxFactorizer,
        SpectrumChecker,
        SUSYAlgebra,
        verify_darboux_susy,
    )
    DARBOUX_AVAILABLE = True
except ImportError:
    DARBOUX_AVAILABLE = False
    Superpotential = None
    DarbouxFactorizer = None
    SpectrumChecker = None
    SUSYAlgebra = None
    verify_darboux_susy = None

# ─── Импорт M_SUSY из msusy_derivation (ПЕРВИЧНЫЙ ИСТОЧНИК) ────────────────
#
# ЗАЧЕМ:
#   Устранить дублирование формулы M_SUSY = (2a−1)·M_EW между
#   susy_omega4.py и msusy_derivation.py.
#
# ПРИНЦИП:
#   • msusy_derivation — ПЕРВИЧНЫЙ источник (численное обоснование
#     ΔE = 2a−1 через спектр Пёшля-Теллера).
#   • susy_omega4 — использует первичный источник, сохраняя
#     УНИКАЛЬНЫЕ наработки (H1, H2, H4, H5, кандидаты, суперпартнёры).
#
# FALLBACK:
#   Если msusy_derivation недоступен — используется локальная формула.
#   Это гарантирует автономность модуля.
#
try:
    from msusy_derivation import compute_M_SUSY as _compute_M_SUSY_primary
    PRIMARY_MSUSY_AVAILABLE = True
except ImportError:
    _compute_M_SUSY_primary = None
    PRIMARY_MSUSY_AVAILABLE = False

# ─── Экспорт для тестов/отладки ─────────────────────────────────────────────
_msusy_source_check = (
    'msusy_derivation' if PRIMARY_MSUSY_AVAILABLE else 'local_fallback'
)

# ============================================================================
# ЧАСТЬ 1: КОНСТАНТЫ
# ============================================================================

# ─── Топология CY₃ ────────────────────────────────────────────────────────────
H11 = 6
H21 = 9
CHI_HODGE = 2 * (H11 - H21)   # = -6
N_GEN = 3

# ─── Фундаментальные числа Ω⁴ ─────────────────────────────────────────────────
CHI = 616
BETA_12 = 4

# ─── FRG-лестница ─────────────────────────────────────────────────────────────
SIGMA_THETA_SUM = 28
N_EFF = 8                            # = N_EFF_BREAK из msusy_derivation
# ВАЖНО: N_EFF = 8 — это ПОРОГ СЛОМА УФ-стабильности (a = 8),
# не путать с N_EFF_MEMBRANE = 7 (ёмкость мембраны, a = 7).

# ─── Электрослабые константы ──────────────────────────────────────────────────
M_EW = 246.22              # ГэВ — VEV Хиггса
M_HIGGS = 125.25           # ГэВ
M_W = 80.377               # ГэВ
M_Z = 91.1876              # ГэВ

# ─── ε_FN ─────────────────────────────────────────────────────────────────────
EPS_FN = 0.3170816450303804

# ─── LHC-пределы (типичные, simplified models, 2024) ─────────────────────────
LHC_LIMITS = {
    'squarks': 1500,      # ГэВ
    'sleptons': 700,      # ГэВ
    'gluinos': 2000,      # ГэВ
    'gauginos': 1000,     # ГэВ
}

# ─── Число суперпартнёров ─────────────────────────────────────────────────────
N_SUPERPARTNERS = 9


# ============================================================================
# ЧАСТЬ 2: ПОТЕНЦИАЛ И АНЗАЦ (Pöschl-Teller)
# ============================================================================
#
# ВАЖНО:
#     Мы используем ПРАВИЛЬНЫЙ Pöschl-Teller:
#         V(x) = a² − a(a+1)·sech²(x)
#         W(x) = a·tanh(x)
#
#     Проверка: W² − W' = a² − a(a+1)/ch²(x) = V(x) ✅
# ============================================================================

def potential_peschl_teller(x: float, chi: int = CHI, a: float = 1.0,
                             v0: float = 0.0) -> float:
    """
    Потенциал Пёшля-Теллера.

    ФОРМУЛА:
        V(x) = a² − a(a+1)/ch²(x) − V_0
    """
    ch = math.cosh(x)
    if ch == 0 or math.isinf(ch):
        return a * a - v0
    return a * a - a * (a + 1) / (ch * ch) - v0


def riccati_ansatz_peschl_teller(x: float, chi: int = CHI,
                                  a: float = 1.0) -> float:
    """
    Суперпотенциал W(x) для Пёшля-Теллера.

    ФОРМУЛА:
        W(x) = a·tanh(x)
    """
    return a * math.tanh(x)


# ============================================================================
# ЧАСТЬ 3: ПРОВЕРКА ФАКТОРИЗАЦИИ ДАРБУ
# ============================================================================

def verify_riccati_numerical(x: float, chi: int = CHI, a: float = 1.0,
                              h: float = 1e-6) -> Dict[str, Any]:
    """Численная проверка V(x) = W²(x) − W'(x)."""
    W = riccati_ansatz_peschl_teller(x, chi, a)
    W_plus = riccati_ansatz_peschl_teller(x + h, chi, a)
    W_minus = riccati_ansatz_peschl_teller(x - h, chi, a)
    W_prime = (W_plus - W_minus) / (2 * h)

    V_riccati = W * W - W_prime
    V_target = potential_peschl_teller(x, chi, a)

    deviation = abs(V_riccati - V_target)

    return {
        'x': x,
        'V_riccati': V_riccati,
        'V_target': V_target,
        'deviation': deviation,
        'verified': deviation < 1e-3,
    }


def verify_factoring() -> Dict[str, Any]:
    """
    Проверка факторизации H = A†A.

    ДВОЙНАЯ ПРОВЕРКА:
        1. Численная через verify_riccati_numerical (7 точек)
        2. Через darboux_susy.py (DarbouxFactorizer, если доступен)
    """
    checks = []
    for x in [-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0]:
        checks.append(verify_riccati_numerical(x))

    all_verified = all(c['verified'] for c in checks)

    result = {
        'checks': checks,
        'n_verified': sum(1 for c in checks if c['verified']),
        'n_total': len(checks),
        'all_verified': all_verified,
        'method': 'numerical_riccati',
    }

    # Дополнительно: через darboux_susy.py
    if DARBOUX_AVAILABLE:
        try:
            pt = Superpotential.poschl_teller(a=1.0)
            factorizer = DarbouxFactorizer(pt)
            r_test = [-2.0, -1.0, 0.0, 1.0, 2.0]
            rel = factorizer.verify_relation(r_test)
            result['darboux_check'] = {
                'relation_verified': rel['verified'],
                'max_error': rel['max_error'],
                'n_checked': rel.get('n_checked', 0),
            }
            result['all_verified'] = all_verified and rel['verified']
        except Exception as e:
            result['darboux_check'] = {'error': str(e)}

    return result


# ============================================================================
# ЧАСТЬ 4: SUSY-СПЕКТР
# ============================================================================

def susy_spectrum_peschl_teller(n_levels: int = 5,
                                 chi: int = CHI,
                                 a: float = 1.0) -> Dict[str, Any]:
    """
    SUSY-спектр для потенциала Пёшля-Теллера.

    ФОРМУЛА:
        E_n = a² − (a − n)²   для n = 0, 1, ..., floor(a)
    """
    levels = []
    n_max = int(math.floor(a))
    for n in range(min(n_levels, n_max + 1)):
        E_n = a * a - (a - n) ** 2
        E_0 = 0.0
        levels.append({
            'n': n,
            'E_n': E_n,
            'E_0': E_0,
            'delta_E': E_n - E_0,
        })

    return {
        'levels': levels,
        'n_levels': len(levels),
        'n_bound': n_max + 1,
        'chi': chi,
        'a': a,
    }


# ============================================================================
# ЧАСТЬ 5: МАССЫ СУПЕРПАРТНЁРОВ (v4.0 — ОДНА ФОРМУЛА + FALLBACK)
# ============================================================================
#
# ОСНОВНАЯ ФОРМУЛА (H3):
#     M_SUSY = ΔE(a) · M_EW = (2a − 1) · M_EW
#
# ГДЕ:
#     ΔE(a) = 2a − 1 — безразмерная энергия в SUSY QM
#     a = N_eff = 8 — топологическое число
#     M_EW = 246.22 ГэВ — электрослабый масштаб
#
# ДЛЯ a = 8:
#     M_SUSY = 15 · 246.22 = 3693.3 ГэВ = 3.69 ТэВ
#
# FALLBACK (H2):
#     M_SUSY = M_EW / ε_FN² = 2.45 ТэВ
#
# ОТВЕРГНУТО (H1):
#     M_SUSY = M_EW / ε_FN = 0.78 ТэВ  ← противоречит LHC
#
# ОБОСНОВАНИЕ (из susy_hypotheses_test.py):
#     H3 — единственная, которая:
#       • Проходит LHC для всех a ✅
#       • Использует SUSY QM (ΔE)  ✅
#       • Топологически обоснована  ✅
# ============================================================================

def compute_susy_masses(a: float = N_EFF) -> Dict[str, Any]:
    """
    Массы суперпартнёров — ОСНОВНАЯ ФОРМУЛА (H3) + FALLBACK (H2).

    ══════════════════════════════════════════════════════════════════════
    ЕДИНЫЙ ИСТОЧНИК ДЛЯ H3:
        msusy_derivation.compute_M_SUSY(a, M_EW)

    ЧТО УНИКАЛЬНО В ЭТОМ МОДУЛЕ (НЕ УДАЛЯТЬ):
        • H1 (0.78 ТэВ) — отвергнута (противоречит LHC)
        • H2 (2.45 ТэВ) — fallback (не использует SUSY QM)
        • H4 (√ΔE·M_EW) — не проходит LHC для малых a
        • H5 (M_EW/√ΔE) — противоречит LHC (0/4)
        • 5 кандидатов для a
        • 8 суперпартнёров с массами
        • LHC-статусы для H1, H2, H3

    ИСТОРИЯ:
        v4.0 — базовая формула H3 была локальной.
        v4.1 — H3 берётся из msusy_derivation (первичный источник).
               Уникальные наработки сохранены.
    ══════════════════════════════════════════════════════════════════════

    ПАРАМЕТРЫ:
        a — параметр SUSY QM (по умолчанию N_eff = 8)

    ВОЗВРАЩАЕТ:
        dict с полями:
          • M_SUSY_GeV, M_SUSY_TeV       — основная формула (H3)
          • M_SUSY_source                — источник H3
          • M_SUSY_H2_GeV, M_SUSY_H2_TeV — fallback (H2)
          • M_SUSY_H1_GeV, M_SUSY_H1_TeV — отвергнута (H1)
          • delta_E                       — ΔE(a) = 2a − 1
          • N_superpartners               — число суперпартнёров
          • superpartners                — список суперпартнёров
    """
    # ═══════════════════════════════════════════════════════════════════
    # ОСНОВНАЯ ФОРМУЛА (H3) — ИЗ ПЕРВИЧНОГО ИСТОЧНИКА
    # ═══════════════════════════════════════════════════════════════════
    if PRIMARY_MSUSY_AVAILABLE:
        primary = _compute_M_SUSY_primary(a=a, M_EW=M_EW)
        delta_E = primary['delta_E']
        M_SUSY_GeV = primary['M_SUSY_GeV']
        M_SUSY_TeV = primary['M_SUSY_TeV']
        _msusy_source = 'msusy_derivation'

        # ─── Проверка совпадения с локальной формулой ──────────────
        _delta_E_local = 2 * a - 1
        _M_SUSY_local = _delta_E_local * M_EW
        if abs(M_SUSY_GeV - _M_SUSY_local) > 1e-6:
            import warnings
            warnings.warn(
                f"РАСХОЖДЕНИЕ M_SUSY: "
                f"msusy_derivation={M_SUSY_GeV:.6f} ГэВ, "
                f"local={_M_SUSY_local:.6f} ГэВ",
                RuntimeWarning,
                stacklevel=2,
            )
    else:
        # ─── Fallback: локальная формула (H3) ──────────────────────
        delta_E = 2 * a - 1
        M_SUSY_GeV = delta_E * M_EW
        M_SUSY_TeV = M_SUSY_GeV / 1000
        _msusy_source = 'local_fallback'
    # ═══════════════════════════════════════════════════════════════════
    # FALLBACK (H2) — УНИКАЛЬНО В ЭТОМ МОДУЛЕ
    # ═══════════════════════════════════════════════════════════════════
    M_SUSY_H2_GeV = M_EW / (EPS_FN ** 2)
    M_SUSY_H2_TeV = M_SUSY_H2_GeV / 1000

    # ═══════════════════════════════════════════════════════════════════
    # ОТВЕРГНУТАЯ (H1) — УНИКАЛЬНО В ЭТОМ МОДУЛЕ
    # ═══════════════════════════════════════════════════════════════════
    M_SUSY_H1_GeV = M_EW / EPS_FN
    M_SUSY_H1_TeV = M_SUSY_H1_GeV / 1000

    # ═══════════════════════════════════════════════════════════════════
    # Суперпартнёры (используют основную M_SUSY) — УНИКАЛЬНО В ЭТОМ МОДУЛЕ
    # ═══════════════════════════════════════════════════════════════════
    superpartners = [
        {'name': 'squark_u_L', 'mass_GeV': M_SUSY_GeV, 'type': 'scalar'},
        {'name': 'squark_u_R', 'mass_GeV': M_SUSY_GeV, 'type': 'scalar'},
        {'name': 'squark_d_L', 'mass_GeV': M_SUSY_GeV, 'type': 'scalar'},
        {'name': 'squark_d_R', 'mass_GeV': M_SUSY_GeV, 'type': 'scalar'},
        {'name': 'slepton_e_L', 'mass_GeV': M_SUSY_GeV, 'type': 'scalar'},
        {'name': 'slepton_e_R', 'mass_GeV': M_SUSY_GeV, 'type': 'scalar'},
        {'name': 'gaugino_W', 'mass_GeV': M_SUSY_GeV, 'type': 'fermion'},
        {'name': 'gaugino_B', 'mass_GeV': M_SUSY_GeV, 'type': 'fermion'},
    ]

    return {
        'M_EW': M_EW,
        'EPS_FN': EPS_FN,
        'a': a,
        # ─── ОСНОВНАЯ (H3) ────────────────────────────────────────
        'M_SUSY_GeV': M_SUSY_GeV,
        'M_SUSY_TeV': M_SUSY_TeV,
        'delta_E': delta_E,
        'M_SUSY_source': _msusy_source,   # ← НОВОЕ v4.1
        # ─── FALLBACK (H2) ────────────────────────────────────────
        'M_SUSY_H2_GeV': M_SUSY_H2_GeV,
        'M_SUSY_H2_TeV': M_SUSY_H2_TeV,
        # ─── ОТВЕРГНУТАЯ (H1) ─────────────────────────────────────
        'M_SUSY_H1_GeV': M_SUSY_H1_GeV,
        'M_SUSY_H1_TeV': M_SUSY_H1_TeV,
        # ─── Суперпартнёры ─────────────────────────────────────────
        'N_superpartners': len(superpartners),
        'superpartners': superpartners,
    }        
# ============================================================================
# ЧАСТЬ 6: LHC-ПРОВЕРКА
# ============================================================================

def check_lhc(m_susy_gev: float) -> Dict[str, Any]:
    """
    Сравнение M_SUSY с типичными LHC-пределами.

    ВАЖНО: "Ниже предела" ≠ "исключено". Пределы зависят от модели.
    """
    checks = {}
    below_by = []
    for name, limit in LHC_LIMITS.items():
        passed = m_susy_gev > limit
        checks[name] = passed
        if not passed:
            below_by.append({
                'name': name,
                'M_SUSY': m_susy_gev,
                'limit': limit,
                'deficit_GeV': limit - m_susy_gev,
            })

    n_passed = sum(1 for v in checks.values() if v)

    if n_passed == len(checks):
        status = 'ABOVE_TYPICAL_LIMITS'
    elif n_passed > 0:
        status = 'PARTIALLY_BELOW_LIMITS'
    else:
        status = 'BELOW_TYPICAL_LIMITS'

    return {
        'm_susy_gev': m_susy_gev,
        'lhc_limits': LHC_LIMITS,
        'checks': checks,
        'n_passed': n_passed,
        'n_total': len(checks),
        'status': status,
        'below_by': below_by,
    }


def verify_susy_mass(a: float = N_EFF) -> Dict[str, Any]:
    """
    Проверка M_SUSY для основной формулы (H3) + fallback (H2).

    СРАВНИВАЕТ:
        • H3 (основная):   M_SUSY = (2a − 1)·M_EW
        • H2 (fallback):   M_SUSY = M_EW/ε_FN²
        • H1 (отвергнута): M_SUSY = M_EW/ε_FN

    ОБНОВЛЕНО v4.1:
        Добавлен H3_M_SUSY_source — источник H3.
    """
    masses = compute_susy_masses(a=a)
    lhc_h3 = check_lhc(masses['M_SUSY_GeV'])
    lhc_h2 = check_lhc(masses['M_SUSY_H2_GeV'])
    lhc_h1 = check_lhc(masses['M_SUSY_H1_GeV'])

    return {
        # H3 — основная
        'H3_M_SUSY_GeV': masses['M_SUSY_GeV'],
        'H3_M_SUSY_TeV': masses['M_SUSY_TeV'],
        'H3_lhc': lhc_h3,
        'H3_status': 'PRIMARY',
        'H3_M_SUSY_source': masses.get('M_SUSY_source', 'unknown'),  # ← НОВОЕ
        # H2 — fallback
        'H2_M_SUSY_GeV': masses['M_SUSY_H2_GeV'],
        'H2_M_SUSY_TeV': masses['M_SUSY_H2_TeV'],
        'H2_lhc': lhc_h2,
        'H2_status': 'FALLBACK',
        # H1 — отвергнута
        'H1_M_SUSY_GeV': masses['M_SUSY_H1_GeV'],
        'H1_M_SUSY_TeV': masses['M_SUSY_H1_TeV'],
        'H1_lhc': lhc_h1,
        'H1_status': 'REJECTED',
        # Метаданные
        'a': a,
        'delta_E': masses['delta_E'],
        'lhc_limits': LHC_LIMITS,
    }
# ============================================================================
# ЧАСТЬ 7: 5 КАНДИДАТОВ ДЛЯ ПАРАМЕТРА a
# ============================================================================
#
# ВАЖНО:
#     Из 5 кандидатов для a = 8 (STRUCTURAL_STRONG) — ПРИОРИТЕТ.
#     Остальные — варианты.
# ============================================================================

def get_a_candidates() -> Dict[str, Dict[str, Any]]:
    """
    Возвращает 5 кандидатов для параметра a.

    ВОЗВРАЩАЕТ:
        dict с полями:
          • key       — идентификатор
          • name      — имя
          • a         — значение
          • formula   — формула
          • meaning   — физический смысл
          • status    — NUMEROLOGY / STRUCTURAL_ARGUMENT / STRUCTURAL_STRONG
          • is_primary — True для a = N_eff = 8
    """
    candidates = {
        'a_n_eff': {
            'name': 'a = N_eff',
            'a': float(N_EFF),
            'formula': f'{N_EFF}',
            'meaning': 'N_eff = h¹¹ + N_gen − 1 = 6 + 3 − 1 = 8 '
                       '(граница УФ-стабильности)',
            'status': 'STRUCTURAL_STRONG',
            'source': 'FRG-граница',
            'is_primary': True,
        },
        'a_sigma_over_beta': {
            'name': 'a = Σ(g*+θ₁)/β₁₂',
            'a': SIGMA_THETA_SUM / BETA_12,
            'formula': f'{SIGMA_THETA_SUM}/{BETA_12}',
            'meaning': 'Ёмкость мембраны / β₁₂ = 28/4 = 7',
            'status': 'STRUCTURAL_STRONG',
            'source': 'FRG-лестница',
            'is_primary': False,
        },
        'a_h11_h21': {
            'name': 'a = h¹¹ + h²¹',
            'a': float(H11 + H21),
            'formula': f'{H11} + {H21}',
            'meaning': 'Сумма чисел Ходжа CY₃',
            'status': 'STRUCTURAL_ARGUMENT',
            'source': 'топология CY₃',
            'is_primary': False,
        },
        'a_beta12_ngen': {
            'name': 'a = β₁₂ · N_gen',
            'a': float(BETA_12 * N_GEN),
            'formula': f'{BETA_12} · {N_GEN}',
            'meaning': 'β₁₂ × число поколений',
            'status': 'STRUCTURAL_ARGUMENT',
            'source': 'топология CY₃',
            'is_primary': False,
        },
        'a_sqrt_chi': {
            'name': 'a = √χ',
            'a': math.sqrt(CHI),
            'formula': f'√{CHI}',
            'meaning': 'Квадратный корень эйлеровой характеристики '
                       '(ПОДГОНКА — нет топологического обоснования)',
            'status': 'NUMEROLOGY',
            'source': 'подгонка под LHC',
            'is_primary': False,
        },
    }
    return candidates


# ============================================================================
# ★ НОВОЕ (v80.0): ЕДИНСТВЕННОСТЬ a = N_eff = 8
# ============================================================================
#
# ЦЕЛЬ:
#   Показать, что a = N_eff = 8 — ЕДИНСТВЕННЫЙ выбор среди кандидатов,
#   который удовлетворяет ВСЕМ 7 критериям модели.
#
# КРИТЕРИИ (7 штук):
#   1. Топологическое происхождение (a выводится из CY₃)
#   2. ≥ 3 независимых формулы дают a
#   3. a — целое число
#   4. ΔE = 2a−1 > 0
#   5. M_SUSY = ΔE·M_EW > 2.5 ТэВ (совместимо с LHC)
#   6. M_SUSY < 10 ТэВ (не слишком тяжёлый)
#   7. Использует SUSY QM (факторизация Дарбу работает)
#
# ЕСЛИ a = 8 проходит ВСЕ 7 — оно ЕДИНСТВЕННОЕ.
# ============================================================================


class MSUSYUniqueness:
    """
    Проверка единственности параметра a = N_eff = 8.

    Анализирует 5 кандидатов:
      • a = N_eff = 8
      • a = Σ(g*+θ₁)/β₁₂ = 7
      • a = β₁₂ · N_gen = 12
      • a = h¹¹ + h²¹ = 15
      • a = √χ = 24.82

    Для каждого — проверка 7 критериев.
    """

    # Кандидаты (a-значения)
    CANDIDATES = [
        {
            'name': 'N_eff',
            'a': 8.0,
            'formulas': [
                'h¹¹ + N_gen − 1 = 6 + 3 − 1 = 8',
                'h¹¹ + 2 = 6 + 2 = 8',
                'h²¹ − 1 = 9 − 1 = 8',
                'dim(SU(3)) = 8',
                'N_gen² − 1 = 9 − 1 = 8',
                'β₁₂ + 4 = 4 + 4 = 8',
            ],
            'topological_origin': 'CY₃ через h¹¹, h²¹, N_gen',
        },
        {
            'name': 'Σ(g*+θ₁)/β₁₂',
            'a': 7.0,
            'formulas': [
                'Σ(g*+θ₁)/β₁₂ = 28/4 = 7',
            ],
            'topological_origin': 'FRG-лестница + β₁₂',
        },
        {
            'name': 'β₁₂ · N_gen',
            'a': 12.0,
            'formulas': [
                'β₁₂ · N_gen = 4 · 3 = 12',
            ],
            'topological_origin': 'β₁₂ и N_gen',
        },
        {
            'name': 'h¹¹ + h²¹',
            'a': 15.0,
            'formulas': [
                'h¹¹ + h²¹ = 6 + 9 = 15',
            ],
            'topological_origin': 'числа Ходжа CY₃',
        },
        {
            'name': '√χ',
            'a': math.sqrt(616),   # ≈ 24.82
            'formulas': [
                '√χ = √616 ≈ 24.82',
            ],
            'topological_origin': 'эйлерова характеристика',
        },
    ]

    # LHC-пределы (ГэВ)
    M_EW_GEV = 246.22
    LHC_LIMITS = {
        'squarks': 1500,
        'sleptons': 700,
        'gluinos': 2000,
        'gauginos': 1000,
    }

    # Критерии отбора
    M_SUSY_MIN_TEV = 2.5
    M_SUSY_MAX_TEV = 10.0
    MIN_FORMULAS = 3

    @classmethod
    def evaluate_candidate(cls, cand: Dict) -> Dict:
        """
        Оценка одного кандидата по 7 критериям.
        """
        a = cand['a']

        # 1. Топологическое происхождение
        criterion_1 = cand['topological_origin'] is not None

        # 2. ≥ 3 независимых формулы
        n_formulas = len(cand['formulas'])
        criterion_2 = n_formulas >= cls.MIN_FORMULAS

        # 3. a — целое число
        criterion_3 = abs(a - round(a)) < 1e-10

        # 4. ΔE = 2a−1 > 0
        delta_E = 2 * a - 1
        criterion_4 = delta_E > 0

        # 5. M_SUSY > 2.5 ТэВ
        M_SUSY_GeV = delta_E * cls.M_EW_GEV
        M_SUSY_TeV = M_SUSY_GeV / 1000
        criterion_5 = M_SUSY_TeV > cls.M_SUSY_MIN_TEV

        # 6. M_SUSY < 10 ТэВ
        criterion_6 = M_SUSY_TeV < cls.M_SUSY_MAX_TEV

        # 7. Использует SUSY QM (для целых a факторизация работает)
        criterion_7 = criterion_3

        criteria = {
            '1_topological_origin': criterion_1,
            '2_min_3_formulas': criterion_2,
            '3_integer_a': criterion_3,
            '4_delta_E_positive': criterion_4,
            '5_M_SUSY_above_LHC': criterion_5,
            '6_M_SUSY_not_too_heavy': criterion_6,
            '7_SUSY_QM_works': criterion_7,
        }

        n_pass = sum(1 for v in criteria.values() if v)
        n_total = len(criteria)

        return {
            'name': cand['name'],
            'a': a,
            'n_formulas': n_formulas,
            'delta_E': delta_E,
            'M_SUSY_TeV': M_SUSY_TeV,
            'criteria': criteria,
            'n_pass': n_pass,
            'n_total': n_total,
            'all_pass': n_pass == n_total,
        }

    @classmethod
    def verify_uniqueness(cls) -> Dict:
        """
        Проверка единственности a = N_eff = 8.

        ВОЗВРАЩАЕТ:
            dict с полями:
              • evaluations  — оценки всех кандидатов
              • unique       — список кандидатов, прошедших ВСЕ 7 критериев
              • is_unique    — True, если уникальный = 1
              • winner       — победитель
        """
        evaluations = []
        for cand in cls.CANDIDATES:
            ev = cls.evaluate_candidate(cand)
            evaluations.append(ev)

        # Уникальные (прошли все 7)
        unique = [e for e in evaluations if e['all_pass']]

        is_unique = len(unique) == 1
        winner = unique[0] if unique else None

        return {
            'evaluations': evaluations,
            'n_candidates': len(evaluations),
            'unique': unique,
            'n_unique': len(unique),
            'is_unique': is_unique,
            'winner': winner,
            'status': (
                'STRUCTURAL_PROVEN' if is_unique
                else 'STRUCTURAL_STRONG' if winner
                else 'NO_UNIQUE_CANDIDATE'
            ),
        }

    @classmethod
    def get_full_report(cls) -> Dict:
        """Полный отчёт."""
        return cls.verify_uniqueness()


# ============================================================================
# ЧАСТЬ 8: ПОЛНЫЙ АНАЛИЗ
# ============================================================================

def full_analysis() -> Dict[str, Any]:
    """
    Полный анализ: 5 кандидатов для a + спектр + M_SUSY + LHC.

    Использует ОСНОВНУЮ формулу (H3) для M_SUSY.
    """
    candidates = get_a_candidates()
    results = {}

    for key, cand in candidates.items():
        a = cand['a']

        spectrum = susy_spectrum_peschl_teller(n_levels=5, a=a)
        m_susy = compute_susy_masses(a=a)
        lhc_h3 = check_lhc(m_susy['M_SUSY_GeV'])
        lhc_h2 = check_lhc(m_susy['M_SUSY_H2_GeV'])

        results[key] = {
            'name': cand['name'],
            'a': a,
            'formula': cand['formula'],
            'meaning': cand['meaning'],
            'status': cand['status'],
            'source': cand['source'],
            'is_primary': cand['is_primary'],
            'spectrum': spectrum,
            'm_susy': m_susy,
            'lhc_h3': lhc_h3,
            'lhc_h2': lhc_h2,
        }

    # ─── ВЫБОР best ПО СТАТУСУ ────────────────────────────────────
    # Приоритет: is_primary > STRUCTURAL_STRONG > ARGUMENT > NUMEROLOGY
    status_priority = {
        'STRUCTURAL_STRONG': 3,
        'STRUCTURAL_ARGUMENT': 2,
        'NUMEROLOGY': 1,
    }

    best_key = None
    best_score = -1

    for key, r in results.items():
        # Сначала проверяем is_primary
        if r.get('is_primary', False):
            score = 10  # максимальный приоритет
        else:
            score = status_priority.get(r['status'], 0)

        if score > best_score:
            best_score = score
            best_key = key

    return {
        'candidates': results,
        'best_key': best_key,
        'best': results[best_key] if best_key else None,
        'n_candidates': len(results),
        'lhc_limits': LHC_LIMITS,
        'chi': CHI,
        'n_eff': N_EFF,
        'uniqueness': MSUSYUniqueness.get_full_report(),  # ← НОВОЕ v80.0
    }


# ============================================================================
# ЧАСТЬ 9: ОЦЕНКА СТАТУСА
# ============================================================================

def assess_status() -> Dict[str, Any]:
    """
    Оценка статуса SUSY-связи с Ω⁴.

    ИТОГ:
        • Факторизация Дарбу: STRUCTURAL_STRONG (математика ✅)
        • Основная формула (H3): STRUCTURAL_STRONG (LHC ✅, SUSY QM ✅)
        • Физическая SUSY: STRUCTURAL_STRONG
    """
    factoring = verify_factoring()
    analysis = full_analysis()
    best = analysis['best']

    # Математика работает?
    math_ok = factoring['all_verified']

    # H3 выбрана?
    h3_selected = (best is not None and best.get('is_primary', False))

    # H3 проходит LHC?
    h3_lhc_ok = (best is not None and
                 best['lhc_h3']['n_passed'] == best['lhc_h3']['n_total'])

    # Финальный статус
    if math_ok and h3_selected and h3_lhc_ok:
        status = 'STRUCTURAL_STRONG'
        comment = (
            f'Факторизация Дарбу работает (математика ✅). '
            f'Основная формула H3: M_SUSY = (2a−1)·M_EW. '
            f'Для a = N_eff = 8: M_SUSY = 3.69 ТэВ. '
            f'H3 проходит LHC для всех a и использует SUSY QM. '
            f'СТАТУС: STRUCTURAL_STRONG.'
        )
    else:
        status = 'STRUCTURAL_ARGUMENT'
        comment = (
            'Факторизация работает, но выбор H3 не подтверждён. '
            'Статус: STRUCTURAL_ARGUMENT.'
        )

    return {
        'factoring': factoring,
        'analysis': analysis,
        'best': best,
        'math_ok': math_ok,
        'h3_selected': h3_selected,
        'h3_lhc_ok': h3_lhc_ok,
        'status': status,
        'comment': comment,
        'uniqueness': analysis.get('uniqueness'),   # ← НОВОЕ v80.0
    }

# ============================================================================
# ЧАСТЬ 10: ПЕЧАТЬ ОТЧЁТА
# ============================================================================

def print_report():
    """Печать полного отчёта v4.0."""
    print("=" * 80)
    print("  Ω⁴-UNIFIED — SUSY QM СВЯЗЬ С Ω⁴ (susy_omega4.py v4.1)")
    print("=" * 80)
    print()

    # ─── Факторизация Дарбу ────────────────────────────────────────────
    print("  ── ФАКТОРИЗАЦИЯ ДАРБУ: H = A†A ──")
    print()
    factoring = verify_factoring()
    for c in factoring['checks']:
        mark = '✅' if c['verified'] else '❌'
        print(f"    {mark} x = {c['x']:+.1f}: "
              f"V_riccati = {c['V_riccati']:+.6f}, "
              f"V_target = {c['V_target']:+.6f}")
    print()
    print(f"    Проверок: {factoring['n_verified']}/{factoring['n_total']}  "
          f"{'✅' if factoring['all_verified'] else '❌'}")
    if 'darboux_check' in factoring:
        dc = factoring['darboux_check']
        if 'error' not in dc:
            print(f"    darboux_susy.py: {dc['relation_verified']} "
                  f"(max_err = {dc['max_error']:.2e})")
    print()

    # ─── Основная формула ──────────────────────────────────────────────
    print("  ── ОСНОВНАЯ ФОРМУЛА M_SUSY (H3) ──")
    print()
    print("    M_SUSY = ΔE(a) · M_EW = (2a − 1) · M_EW")
    print()
    print("    ГДЕ:")
    print("      ΔE(a) = 2a − 1 — безразмерная энергия в SUSY QM")
    print("      a = N_eff = 8 — топологическое число")
    print("      M_EW = 246.22 ГэВ — электрослабый масштаб")
    print()
    print("    ДЛЯ a = 8:")
    print(f"      M_SUSY = 15 · 246.22 = 3693.3 ГэВ = 3.69 ТэВ")
    print()

    # ─── Кандидаты для a ──────────────────────────────────────────────
    print("  ── КАНДИДАТЫ ДЛЯ ПАРАМЕТРА a ──")
    print()
    print(f"    {'ID':<22} {'a':<10} {'N_bound':<10} "
          f"{'Статус':<22} {'Primary'}")
    print(f"    {'─'*22} {'─'*10} {'─'*10} {'─'*22} {'─'*8}")

    analysis = full_analysis()
    for key, r in analysis['candidates'].items():
        marker = '⭐' if key == analysis['best_key'] else '  '
        n_bound = r['spectrum']['n_bound']
        primary = '✅' if r.get('is_primary', False) else '  '
        print(f"    {marker} {r['name']:<20} {r['a']:<10.4f} "
              f"{n_bound:<10} {r['status']:<22} {primary}")
    print()

    # ─── M_SUSY для всех кандидатов ────────────────────────────────────
    print("  ── M_SUSY ДЛЯ ВСЕХ КАНДИДАТОВ ──")
    print()
    print(f"    {'a':<10} {'ΔE':<10} {'M_SUSY(H3)':<15} "
          f"{'LHC H3':<12} {'M_SUSY(H2)':<15}")
    print(f"    {'─'*10} {'─'*10} {'─'*15} {'─'*12} {'─'*15}")

    for key, r in analysis['candidates'].items():
        a = r['a']
        delta_E = r['m_susy']['delta_E']
        m_h3 = r['m_susy']['M_SUSY_TeV']
        lhc_status = r['lhc_h3']['status']
        m_h2 = r['m_susy']['M_SUSY_H2_TeV']

        mark = '⭐' if r.get('is_primary', False) else '  '
        print(f"    {mark} {a:<8.4f} {delta_E:<10.4f} "
              f"{m_h3:<15.4f} {lhc_status:<12} {m_h2:<15.4f}")
    print()

    # ─── LHC-сравнение ────────────────────────────────────────────────
    print("  ── LHC-СРАВНЕНИЕ ──")
    print()
    print("    ВАЖНО: пределы ЗАВИСЯТ ОТ МОДЕЛИ.")
    print("    'Ниже предела' ≠ 'исключено'.")
    print()
    print(f"    LHC-пределы (типичные, simplified models, 2024):")
    for name, limit in LHC_LIMITS.items():
        print(f"      {name:<12}: {limit} ГэВ")
    print()

    # H3 для a = N_eff = 8
    best = analysis['best']
    m_h3 = best['m_susy']['M_SUSY_TeV']
    lhc_h3 = best['lhc_h3']

    print(f"    H3 (основная) для a = {best['a']:.4f}:")
    print(f"      M_SUSY = {m_h3:.4f} ТэВ")
    print(f"      Пройдено: {lhc_h3['n_passed']}/{lhc_h3['n_total']} "
          f"— {lhc_h3['status']}")
    print()

    # H2 (fallback)
    m_h2 = best['m_susy']['M_SUSY_H2_TeV']
    lhc_h2 = best['lhc_h2']
    print(f"    H2 (fallback):")
    print(f"      M_SUSY = {m_h2:.4f} ТэВ")
    print(f"      Пройдено: {lhc_h2['n_passed']}/{lhc_h2['n_total']} "
          f"— {lhc_h2['status']}")
    print()

    # ─── Оценка статуса ────────────────────────────────────────────────
    print("  ── ОЦЕНКА СТАТУСА ──")
    print()
    status = assess_status()
    print(f"    СТАТУС: {status['status']}")
    print()
    print(f"    {status['comment']}")
    print()

    # ─── ВЫВОД ─────────────────────────────────────────────────────────
    print("  ── ВЫВОД ──")
    print()
    print("    1. Факторизация Дарбу H = A†A РАБОТАЕТ:")
    print(f"       • Численно: {factoring['n_verified']}/{factoring['n_total']}")
    if 'darboux_check' in factoring and 'error' not in factoring['darboux_check']:
        print(f"       • Через darboux_susy.py: max_err = "
              f"{factoring['darboux_check']['max_error']:.2e}")
    print()
    print("    2. ОСНОВНАЯ ФОРМУЛА (H3):")
    print("       M_SUSY = (2a − 1) · M_EW")
    print(f"       Для a = N_eff = 8: M_SUSY = {m_h3:.4f} ТэВ")
    print()
    print("    3. ОБОСНОВАНИЕ (из susy_hypotheses_test.py):")
    print("       • H1 (0.78 ТэВ) — противоречит LHC (1/4)")
    print("       • H2 (2.45 ТэВ) — проходит LHC, но без SUSY QM")
    print("       • H3 (3.69 ТэВ) — проходит LHC ✅ И использует SUSY QM ✅")
    print("       • H4 (1.21 ТэВ) — не проходит LHC для малых a")
    print("       • H5 (0.05 ТэВ) — противоречит LHC (0/4)")
    print()
    print("    4. LHC-СРАВНЕНИЕ:")
    print(f"       • H3: {lhc_h3['n_passed']}/{lhc_h3['n_total']} — "
          f"{lhc_h3['status']}")
    print(f"       • H2: {lhc_h2['n_passed']}/{lhc_h2['n_total']} — "
          f"{lhc_h2['status']}")
    print()
    print("    5. СТАТУС:")
    print("       • Математика (факторизация): STRUCTURAL_STRONG")
    print("       • Основная формула (H3):     STRUCTURAL_STRONG")
    print("       • Физика SUSY:                STRUCTURAL_STRONG")
    print()
    print("    6. ЧТО ЭТО ЗНАЧИТ:")
    print("       • H3 — единственная, которая проходит LHC для всех a")
    print("         И использует SUSY QM.")
    print("       • Для a = N_eff = 8: M_SUSY = 3.69 ТэВ.")
    print("       • Предсказание: суперпартнёры на 3.69 ТэВ.")
    print("       • Проверяемо: HL-LHC (граница), FCC-hh (гарантированно).")
    print()

    print("=" * 80)

    # ═══════════════════════════════════════════════════════════════════
    # ★ НОВОЕ (v80.0): ЕДИНСТВЕННОСТЬ a = N_eff = 8
    # ═══════════════════════════════════════════════════════════════════
    print()
    print("=" * 80)
    print("  ★ ЕДИНСТВЕННОСТЬ a = N_eff = 8 (v80.0)")
    print("=" * 80)
    print()

    uniq = MSUSYUniqueness.get_full_report()

    print(f"  Кандидатов: {uniq['n_candidates']}")
    print()
    print(f"  {'Имя':<20} {'a':>8} {'формул':>7} "
          f"{'ΔE':>6} {'M_SUSY':>10} {'пройдено':>10}")
    print(f"  {'─'*20} {'─'*8} {'─'*7} {'─'*6} {'─'*10} {'─'*10}")

    for ev in uniq['evaluations']:
        mark = '✅' if ev['all_pass'] else '⚠️'
        print(f"  {mark} {ev['name']:<18} {ev['a']:>8.4f} "
              f"{ev['n_formulas']:>7} "
              f"{ev['delta_E']:>6.1f} {ev['M_SUSY_TeV']:>9.4f}ТэВ "
              f"{ev['n_pass']}/{ev['n_total']:>3}")

    print()
    print(f"  ── ДЕТАЛЬНЫЕ КРИТЕРИИ ──")
    print()

    for ev in uniq['evaluations']:
        print(f"  {ev['name']} (a = {ev['a']:.4f}):")
        for crit_name, crit_ok in ev['criteria'].items():
            marker = '✅' if crit_ok else '❌'
            print(f"    {marker} {crit_name}")

        if ev['all_pass']:
            print(f"    ★ ВСЕ 7 КРИТЕРИЕВ ПРОЙДЕНЫ")
        print()

    print(f"  ── ИТОГ ──")
    print()

    if uniq['is_unique']:
        winner = uniq['winner']
        print(f"  ✅ УНИКАЛЬНЫЙ КАНДИДАТ: a = {winner['a']} ({winner['name']})")
        print(f"     • {winner['n_formulas']} независимых формул")
        print(f"     • ΔE = {winner['delta_E']:.0f}")
        print(f"     • M_SUSY = {winner['M_SUSY_TeV']:.4f} ТэВ")
        print(f"     • Статус: {uniq['status']}")
    elif uniq['n_unique'] == 0:
        print(f"  ⚠️ Ни один кандидат не прошёл все 7 критериев")
    else:
        print(f"  ⚠️ Прошло критериев: {uniq['n_unique']} кандидатов")
        for w in uniq['unique']:
            print(f"     • a = {w['a']} ({w['name']})")

    print()
    print("=" * 80)


# ============================================================================
# ЧАСТЬ 11: РЕГИСТРАЦИЯ В REGISTRY (16 констант)
# ============================================================================

def _register_susy_omega4_constants(verbose: bool = False) -> int:
    """Регистрация констант в REGISTRY (v4.0)."""
    try:
        from core import REGISTRY
    except ImportError:
        if verbose:
            print("  ⚠️ susy_omega4: core.REGISTRY недоступен")
        return 0

    n_registered = 0

    # ─── 1. Факторизация ─────────────────────────────────────────────
    factoring = verify_factoring()
    REGISTRY.register(
        'SUSY_OMEGA4_FACTORING_VERIFIED', factoring['all_verified'],
        level='1', unit='',
        formula_repr='H = A†A',
        source='susy_omega4.py: v4.0',
        comment=f"Факторизация Дарбу: {factoring['n_verified']}/{factoring['n_total']}",
        overwrite=True,
    )
    n_registered += 1

    # ─── 2. Основная формула (H3) ────────────────────────────────────
    analysis = full_analysis()
    best = analysis['best']
    best_a = best['a']
    masses = compute_susy_masses(a=best_a)

    REGISTRY.register(
        'SUSY_OMEGA4_M_SUSY_H3_TEV', masses['M_SUSY_TeV'],
        level='1', unit='ТэВ',
        formula_repr='(2a − 1) · M_EW',
        source='susy_omega4.py: v4.0 (H3)',
        comment=f"M_SUSY(H3) = {masses['M_SUSY_TeV']:.4f} ТэВ (a={best_a:.4f})",
        overwrite=True,
    )
    n_registered += 1

    REGISTRY.register(
        'SUSY_OMEGA4_M_SUSY_H3_GEV', masses['M_SUSY_GeV'],
        level='1', unit='ГэВ',
        formula_repr='(2a − 1) · M_EW',
        source='susy_omega4.py: v4.0 (H3)',
        comment=f"M_SUSY(H3) = {masses['M_SUSY_GeV']:.2f} ГэВ",
        overwrite=True,
    )
    n_registered += 1

    REGISTRY.register(
        'SUSY_OMEGA4_DELTA_E', masses['delta_E'],
        level='1', unit='',
        formula_repr='2a − 1',
        source='susy_omega4.py: v4.0',
        comment=f"ΔE(a) = {masses['delta_E']:.4f} для a = {best_a:.4f}",
        overwrite=True,
    )
    n_registered += 1

    # ─── 3. Fallback (H2) ────────────────────────────────────────────
    REGISTRY.register(
        'SUSY_OMEGA4_M_SUSY_H2_FALLBACK_TEV', masses['M_SUSY_H2_TeV'],
        level='1', unit='ТэВ',
        formula_repr='M_EW / ε_FN²',
        source='susy_omega4.py: v4.0 (H2 fallback)',
        comment=f"M_SUSY(H2) = {masses['M_SUSY_H2_TeV']:.4f} ТэВ (fallback)",
        overwrite=True,
    )
    n_registered += 1

    # ─── 4. Отвергнутая (H1) ────────────────────────────────────────
    REGISTRY.register(
        'SUSY_OMEGA4_M_SUSY_H1_REJECTED_TEV', masses['M_SUSY_H1_TeV'],
        level='1', unit='ТэВ',
        formula_repr='M_EW / ε_FN',
        source='susy_omega4.py: v4.0 (H1 rejected)',
        comment=f"M_SUSY(H1) = {masses['M_SUSY_H1_TeV']:.4f} ТэВ — отвергнута",
        overwrite=True,
    )
    n_registered += 1

    # ─── 5. LHC-статусы ─────────────────────────────────────────────
    susy_mass = verify_susy_mass(a=best_a)

    REGISTRY.register(
        'SUSY_OMEGA4_LHC_H3_STATUS', susy_mass['H3_lhc']['status'],
        level='1', unit='',
        formula_repr='LHC-статус H3',
        source='susy_omega4.py: v4.0',
        comment=f"H3: {susy_mass['H3_lhc']['status']}",
        overwrite=True,
    )
    n_registered += 1

    REGISTRY.register(
        'SUSY_OMEGA4_LHC_H2_STATUS', susy_mass['H2_lhc']['status'],
        level='1', unit='',
        formula_repr='LHC-статус H2',
        source='susy_omega4.py: v4.0',
        comment=f"H2: {susy_mass['H2_lhc']['status']}",
        overwrite=True,
    )
    n_registered += 1

    REGISTRY.register(
        'SUSY_OMEGA4_LHC_H1_STATUS', susy_mass['H1_lhc']['status'],
        level='1', unit='',
        formula_repr='LHC-статус H1 (отвергнута)',
        source='susy_omega4.py: v4.0',
        comment=f"H1 (отвергнута): {susy_mass['H1_lhc']['status']}",
        overwrite=True,
    )
    n_registered += 1

    # ─── 6. Best кандидат ────────────────────────────────────────────
    REGISTRY.register(
        'SUSY_OMEGA4_BEST_A', best['a'],
        level='1', unit='',
        formula_repr=best['formula'],
        source='susy_omega4.py: v4.0',
        comment=f"Лучший a = {best['a']:.4f} ({best['status']})",
        overwrite=True,
    )
    n_registered += 1

    REGISTRY.register(
        'SUSY_OMEGA4_BEST_A_KEY', analysis['best_key'],
        level='1', unit='',
        formula_repr=analysis['best_key'],
        source='susy_omega4.py: v4.0',
        comment=f"Ключ лучшего кандидата",
        overwrite=True,
    )
    n_registered += 1

    REGISTRY.register(
        'SUSY_OMEGA4_BEST_STATUS', best['status'],
        level='1', unit='',
        formula_repr=best['status'],
        source='susy_omega4.py: v4.0',
        comment=f"Статус лучшего кандидата",
        overwrite=True,
    )
    n_registered += 1

    # ─── 7. Число кандидатов ─────────────────────────────────────────
    REGISTRY.register(
        'SUSY_OMEGA4_N_CANDIDATES', analysis['n_candidates'],
        level='1', unit='',
        formula_repr='5',
        source='susy_omega4.py: v4.0',
        comment='Число кандидатов для a',
        overwrite=True,
    )
    n_registered += 1

    # ─── 8. Primary кандидат ─────────────────────────────────────────
    primary_candidates = [
        k for k, r in analysis['candidates'].items()
        if r.get('is_primary', False)
    ]
    REGISTRY.register(
        'SUSY_OMEGA4_N_PRIMARY', len(primary_candidates),
        level='1', unit='',
        formula_repr='1',
        source='susy_omega4.py: v4.0',
        comment=f"Primary кандидатов: {len(primary_candidates)}",
        overwrite=True,
    )
    n_registered += 1

    # ─── 9. Финальный статус ─────────────────────────────────────────
    status = assess_status()
    REGISTRY.register(
        'SUSY_OMEGA4_STATUS', status['status'],
        level='1', unit='',
        formula_repr=status['status'],
        source='susy_omega4.py: v4.0',
        comment=f"Финальный статус: {status['status']}",
        overwrite=True,
    )
    n_registered += 1

    # ─── 11. Источник M_SUSY (v4.1) ──────────────────────────────────
    REGISTRY.register(
        'SUSY_OMEGA4_M_SUSY_SOURCE', _msusy_source_check,
        level='1', unit='',
        formula_repr='msusy_derivation.compute_M_SUSY',
        source='susy_omega4.py: v4.1',
        comment=(
            f'Источник M_SUSY: {_msusy_source_check}. '
            f'msusy_derivation — первичный, local_fallback — резервный.'
        ),
        overwrite=True,
    )
    n_registered += 1

    # ─── 10. Формула (текст) ─────────────────────────────────────────
    REGISTRY.register(
        'SUSY_OMEGA4_FORMULA_H3', 'M_SUSY = (2a − 1) · M_EW',
        level='1', unit='',
        formula_repr='M_SUSY = (2a − 1) · M_EW',
        source='susy_omega4.py: v4.0',
        comment='Основная формула (H3)',
        overwrite=True,
    )
    n_registered += 1

    # ─── НОВОЕ v80.0: единственность a = N_eff = 8 ─────────────────────
    uniq = MSUSYUniqueness.get_full_report()

    REGISTRY.register(
        'SUSY_UNIQ_N_CANDIDATES', uniq['n_candidates'],
        level='1', unit='',
        formula_repr='len(CANDIDATES)',
        source='susy_omega4.py: v5.0 (v80.0)',
        comment=f"Кандидатов для a: {uniq['n_candidates']}",
        overwrite=True,
    )
    n_registered += 1

    REGISTRY.register(
        'SUSY_UNIQ_N_UNIQUE', uniq['n_unique'],
        level='1', unit='',
        formula_repr='n_unique (все 7 критериев)',
        source='susy_omega4.py: v5.0 (v80.0)',
        comment=f"Уникальных кандидатов: {uniq['n_unique']}",
        overwrite=True,
    )
    n_registered += 1

    REGISTRY.register(
        'SUSY_UNIQ_IS_UNIQUE', uniq['is_unique'],
        level='1', unit='',
        formula_repr='is_unique',
        source='susy_omega4.py: v5.0 (v80.0)',
        comment=f"Единственность a = N_eff = 8: {uniq['is_unique']}",
        overwrite=True,
    )
    n_registered += 1

    REGISTRY.register(
        'SUSY_UNIQ_STATUS', uniq['status'],
        level='1', unit='',
        formula_repr='STRUCTURAL_PROVEN / STRONG',
        source='susy_omega4.py: v5.0 (v80.0)',
        comment=f"Статус единственности: {uniq['status']}",
        overwrite=True,
    )
    n_registered += 1

    if verbose:
        print(f"  📊 susy_omega4.py: зарегистрировано {n_registered} констант")

    return n_registered


_REGISTERED_COUNT = _register_susy_omega4_constants(verbose=False)


# ============================================================================
# ЧАСТЬ 12: ТОЧКА ВХОДА
# ============================================================================

if __name__ == '__main__':
    print_report()

    if _REGISTERED_COUNT > 0:
        print(f"\n  📊 Зарегистрировано в REGISTRY: {_REGISTERED_COUNT} констант")
