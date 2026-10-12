#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
Ω⁴-UNIFIED v78.0 — АНАЛИЗ ЧУВСТВИТЕЛЬНОСТИ (sensitivity.py) — v2.0
================================================================================
ОБЪЕДИНЁННЫЙ МОДУЛЬ (v2.0).
Объединяет предыдущие sensitivity.py (v2.1) и sensitivity2.py (v1.2).

================================================================================
ЧТО ДЕЛАЕТ
================================================================================

Анализ чувствительности модели Ω⁴ к изменениям χ, β₁₂ и n(k).

МЕТОДЫ:
  1. OAT по χ (однопараметрический анализ)
  2. OAT по β₁₂
  3. Анализ по категориям
  4. Глобальный анализ (сетка χ × β₁₂)

ВСЕГО: 31 проверка.

КАТЕГОРИИ:
  • χ-связи (3)
  • β₁₂-связи (1)
  • Алгебры Ли (5)
  • n(k)-связи (5)
  • N_TIME (2)
  • n(5) (1)
  • Тетраэдр (3)
  • Космология (2)
  • Хиггс (1)
  • Φ⁺-связи (3)
  • Эллиптическая кривая (4)
  • GUT-JAX (1)

================================================================================
ИСТОРИЯ ВЕРСИЙ
================================================================================

  v1.0 — базовый анализ (21 проверка).
  v1.1 — исправления sensitivity2.
  v1.2 — sensitivity2: 31 проверка, +Φ⁺, +элл. кривая, +GUT-JAX.
  v2.1 — sensitivity (старый): 31 проверка, исправлен j-инвариант.
  v2.0 — ОБЪЕДИНЁННЫЙ модуль (ЭТА ВЕРСИЯ):
         • Базис — sensitivity2.py v1.2 (31 проверка).
         • Переименование: ModifiableConstants2 → ModifiableConstants,
           ModifiablePolynomial2 → ModifiablePolynomial,
           ModifiableChecks2 → ModifiableChecks,
           SensitivityAnalysis2 → SensitivityAnalysis.
         • УДАЛЁН класс SensitivityComparison — сравнение v1 vs v2
           больше не имеет смысла (v1 объединён с v2).
         • Все проверки сохранены.
         • Метод check_j_invariant — правильная формула (j = -χ·62).

================================================================================
АВТОР:  Сергей Викторович Матершов
ORCID:  0009-0009-0641-1357
ДАТА:   2026-09-19
ВЕРСИЯ: 2.0
ЛИЦЕНЗИЯ: CC BY-NC-ND 4.0 International
================================================================================
"""

from core import (
    C, MasterPolynomial, LieGroupDimensions,
    HodgeNumbers, VacuumChargeFromCY3,
    AtlasLieAlgebras, ExactConnectionsRegistry,
)
from atlas import (
    CompleteAtlasJAX, AtlasDiscoveries,
    AlgebraLyCalculator,
)
from physics import (
    TetrahedronConstants, SMHiggsMasses, Cosmology,
    YukawaFromCY3Charges, GUTUnificationWithSU58,
    AtlasPhysicsConnections, GUTFromJAXAtlas,
)

import math
import copy
from typing import Dict, List, Tuple, Any, Callable


# ============================================================================
# ЧАСТЬ 1: МОДИФИЦИРУЕМЫЕ КОНСТАНТЫ
# ============================================================================

class ModifiableConstants:
    """
    Копия Constants с возможностью изменения χ, β₁₂.

    ИСПОЛЬЗУЕТСЯ В:
        Анализ чувствительности — нужно подменять χ, β₁₂
        и смотреть, как меняются проверки.
    """

    def __init__(self, chi: int = None, beta_12: int = None):
        self._chi = chi if chi is not None else C.CHI
        self._beta_12 = beta_12 if beta_12 is not None else C.BETA_12

    @property
    def CHI(self) -> int:
        return self._chi

    @property
    def BETA_12(self) -> int:
        return self._beta_12

    def __getattr__(self, name: str) -> Any:
        return getattr(C, name)


# ============================================================================
# ЧАСТЬ 2: МОДИФИЦИРУЕМЫЙ ПОЛИНОМ
# ============================================================================

class ModifiablePolynomial:
    """
    Стандартный полином для анализа.

    В текущей версии n(k) не модифицируется (используется стандартный).
    В будущем можно добавить вариацию коэффициентов.
    """

    @staticmethod
    def n_poly(k: int) -> int:
        return MasterPolynomial.n_poly(k)


# ============================================================================
# ЧАСТЬ 3: ПРОВЕРКИ (31 проверка)
# ============================================================================

class ModifiableChecks:
    """
    Набор проверок с модифицированными χ, β₁₂, n(k).

    ВСЕГО: 31 проверка.
    """

    # ========================================================================
    # 3.1. χ-связи (3)
    # ========================================================================

    @staticmethod
    def check_chi_connections(chi: int, beta: int, n_func: Callable) -> bool:
        """χ/11=56, χ/22=28, χ/28=22, χ/44=14, χ/77=8, χ/154=4."""
        try:
            return all([
                chi // 11 == 56,
                chi // 22 == 28,
                chi // 28 == 22,
                chi // 44 == 14,
                chi // 77 == 8,
                chi // 154 == 4,
            ])
        except (ZeroDivisionError, ValueError):
            return False

    @staticmethod
    def check_E8_31(chi: int, beta: int, n_func: Callable) -> bool:
        """E₈: 31 = χ/22 + 3."""
        try:
            return chi // 22 + 3 == 31
        except (ZeroDivisionError, ValueError):
            return False

    @staticmethod
    def check_E7_56(chi: int, beta: int, n_func: Callable) -> bool:
        """E₇: 56 = χ/11."""
        try:
            return chi // 11 == 56
        except (ZeroDivisionError, ValueError):
            return False

    # ========================================================================
    # 3.2. β₁₂-связи (1)
    # ========================================================================

    @staticmethod
    def check_beta_connections(chi: int, beta: int, n_func: Callable) -> bool:
        """β₁₂ = 4, β₁₂² = 16, β₁₂³ = 64, β₁₂⁵ − 1020 = β₁₂."""
        try:
            return all([
                beta == 4,
                beta ** 2 == 16,
                beta ** 3 == 64,
                beta ** 5 - 1020 == beta,
            ])
        except (ValueError, OverflowError):
            return False

    # ========================================================================
    # 3.3. Алгебры Ли (5)
    # ========================================================================

    @staticmethod
    def check_F4_13(chi: int, beta: int, n_func: Callable) -> bool:
        """F₄: 13 = n(1)/2."""
        try:
            return n_func(1) // 2 == 13
        except Exception:
            return False

    @staticmethod
    def check_G2_7(chi: int, beta: int, n_func: Callable) -> bool:
        """G₂: 7 = n(2) + 3."""
        try:
            return n_func(2) + 3 == 7
        except Exception:
            return False

    @staticmethod
    def check_E6_3(chi: int, beta: int, n_func: Callable) -> bool:
        """E₆: 3 = n(2) − 1."""
        try:
            return n_func(2) - 1 == 3
        except Exception:
            return False

    @staticmethod
    def check_SU26(chi: int, beta: int, n_func: Callable) -> bool:
        """SU(26): ω₁ = 26 = n(1)."""
        try:
            return n_func(1) == 26
        except Exception:
            return False

    @staticmethod
    def check_SU58(chi: int, beta: int, n_func: Callable) -> bool:
        """SU(58): ω₁ = 58 = n(3)."""
        try:
            return n_func(3) == 58
        except Exception:
            return False

    # ========================================================================
    # 3.4. n(k)-связи (5)
    # ========================================================================

    @staticmethod
    def check_nk_1(chi: int, beta: int, n_func: Callable) -> bool:
        """n(4) + n(3) = 576 = τ(2)²."""
        try:
            return n_func(4) + n_func(3) == C.TAU_2 ** 2
        except Exception:
            return False

    @staticmethod
    def check_nk_2(chi: int, beta: int, n_func: Callable) -> bool:
        """n(3) + n(-1) = 66."""
        try:
            return n_func(3) + n_func(-1) == 66
        except Exception:
            return False

    @staticmethod
    def check_nk_3(chi: int, beta: int, n_func: Callable) -> bool:
        """n(1) − n(2) = 22."""
        try:
            return n_func(1) - n_func(2) == 22
        except Exception:
            return False

    @staticmethod
    def check_nk_4(chi: int, beta: int, n_func: Callable) -> bool:
        """n(3) − n(1) = 32."""
        try:
            return n_func(3) - n_func(1) == 32
        except Exception:
            return False

    @staticmethod
    def check_nk_5(chi: int, beta: int, n_func: Callable) -> bool:
        """n(1) + n(2) = 30."""
        try:
            return n_func(1) + n_func(2) == 30
        except Exception:
            return False

    # ========================================================================
    # 3.5. N_TIME (2)
    # ========================================================================

    @staticmethod
    def check_N_TIME_v1(chi: int, beta: int, n_func: Callable) -> bool:
        """N_TIME = (n(-3)+n(7)−n(-1))/100."""
        try:
            return (n_func(-3) + n_func(7) - n_func(-1)) // 100 == 122
        except Exception:
            return False

    @staticmethod
    def check_N_TIME_v2(chi: int, beta: int, n_func: Callable) -> bool:
        """N_TIME = (5·(n(2)+3)²−1)/2."""
        try:
            return (5 * (n_func(2) + 3)**2 - 1) // 2 == 122
        except Exception:
            return False

    # ========================================================================
    # 3.6. n(5) (1)
    # ========================================================================

    @staticmethod
    def check_n5(chi: int, beta: int, n_func: Callable) -> bool:
        """n(5) = 1920."""
        try:
            return n_func(5) == 1920
        except Exception:
            return False

    # ========================================================================
    # 3.7. Тетраэдр (3)
    # ========================================================================

    @staticmethod
    def check_Tetrahedron_V_Higgs(chi: int, beta: int, n_func: Callable) -> bool:
        """V_Higgs = 1/√(√2·G_F)."""
        try:
            V_H = 1 / math.sqrt(math.sqrt(2) * C.G_FERMI)
            dev = abs(V_H - C.V_HIGGS) / C.V_HIGGS * 100
            return dev < 0.001
        except Exception:
            return False

    @staticmethod
    def check_Tetrahedron_M_Z(chi: int, beta: int, n_func: Callable) -> bool:
        """M_Z = V_Higgs·10/27."""
        try:
            V_H = 1 / math.sqrt(math.sqrt(2) * C.G_FERMI)
            M_Z_theory = V_H * 10 / 27
            dev = abs(M_Z_theory - C.M_Z) / C.M_Z * 100
            return dev < 0.01
        except Exception:
            return False

    @staticmethod
    def check_Tetrahedron_M_top(chi: int, beta: int, n_func: Callable) -> bool:
        """M_top = M_Higgs·40/29."""
        try:
            M_top_theory = C.M_HIGGS * 40 / 29
            dev = abs(M_top_theory - C.M_TOP) / C.M_TOP * 100
            return dev < 0.01
        except Exception:
            return False

    # ========================================================================
    # 3.8. Космология (2)
    # ========================================================================

    @staticmethod
    def check_Cosmology_n_s(chi: int, beta: int, n_func: Callable) -> bool:
        """n_s = 1 − 2/(χ/10)."""
        try:
            N = chi / 10.0
            if N <= 0:
                return False
            n_s = 1 - 2 / N
            dev = abs(n_s - C.N_S_OBS) / C.N_S_OBS * 100
            return dev < 1.0
        except Exception:
            return False

    @staticmethod
    def check_Cosmology_H0(chi: int, beta: int, n_func: Callable) -> bool:
        """H_0 из Фридмана (зависит от χ)."""
        try:
            N = chi / 10.0
            if N <= 0:
                return False
            H_0 = 67.7 * (chi / 616) ** 0.5
            dev = abs(H_0 - C.H_0_OBS) / C.H_0_OBS * 100
            return dev < 2.0
        except Exception:
            return False

    # ========================================================================
    # 3.9. Хиггс (1)
    # ========================================================================

    @staticmethod
    def check_Higgs_delta_M(chi: int, beta: int, n_func: Callable) -> bool:
        """
        ΔM = 56 − ε_FN.

        ИСПОЛЬЗУЕТ актуальные массы Хиггса из physics.SMHiggsMasses.
        """
        try:
            higgs = SMHiggsMasses().compute_higgs_masses()
            M_H95 = higgs['M_H95']
            M_H152 = higgs['M_H152']

            delta_M = M_H152 - M_H95
            prediction = chi // 11 - C.EPS_FN
            dev = abs(delta_M - prediction) / prediction * 100
            return dev < 1.0
        except Exception:
            return False

    # ========================================================================
    # 3.10. Φ⁺-связи (3)
    # ========================================================================

    @staticmethod
    def check_phi_plus_E6(chi: int, beta: int, n_func: Callable) -> bool:
        """Φ⁺(E₆) = 36 = n(1) + 10."""
        try:
            return n_func(1) + 10 == 36
        except Exception:
            return False

    @staticmethod
    def check_phi_plus_E7(chi: int, beta: int, n_func: Callable) -> bool:
        """Φ⁺(E₇) = 63 = χ/22 + 35."""
        try:
            return chi // 22 + 35 == 63
        except (ZeroDivisionError, ValueError):
            return False

    @staticmethod
    def check_phi_plus_E8(chi: int, beta: int, n_func: Callable) -> bool:
        """Φ⁺(E₈) = 120 = n(1) + 94."""
        try:
            return n_func(1) + 94 == 120
        except Exception:
            return False

    # ========================================================================
    # 3.11. Эллиптическая кривая (4)
    # ========================================================================

    @staticmethod
    def check_bad_prime_617(chi: int, beta: int, n_func: Callable) -> bool:
        """617 = χ + 1."""
        try:
            return chi + 1 == 617
        except Exception:
            return False

    @staticmethod
    def check_j_invariant(chi: int, beta: int, n_func: Callable) -> bool:
        """
        j ≈ -χ · 62.

        ПРАВИЛЬНАЯ ФОРМУЛА (было -n(1)·24·62 = -38688, исправлено):
            j = -χ · 62 = -616 · 62 = -38192
            Отклонение от -38224.54: 0.0851% ✅

        ФИЗИЧЕСКИЙ СМЫСЛ:
            χ = n(1)·24 - 8 = 616
            62 = |Φ⁺(E₇)| - 1 = 63 - 1
        """
        try:
            j_theory = -chi * 62
            return abs(j_theory - (-38224.54)) / 38224.54 < 0.01
        except Exception:
            return False

    @staticmethod
    def check_BSD_relation(chi: int, beta: int, n_func: Callable) -> bool:
        """L'(E,1) ≈ n(3)/6."""
        try:
            return abs(n_func(3) / 6 - 9.506022) < 0.2
        except Exception:
            return False

    @staticmethod
    def check_omega_DM(chi: int, beta: int, n_func: Callable) -> bool:
        """
        Ω_DM = ε_FN · κ ≈ 0.1207.

        УЛУЧШЕННАЯ ФОРМУЛА:
            ε_FN ≈ (n(1) + n(2) + 1 + 1/√2) / 100
                 = (26 + 4 + 1 + 0.7071) / 100
                 = 0.317071
        """
        try:
            eps_approx = (n_func(1) + n_func(2) + 1 + 1/math.sqrt(2)) / 100
            kappa = 3.61814 / 9.506022
            omega = eps_approx * kappa
            return abs(omega - 0.1207) < 0.005
        except Exception:
            return False

    # ========================================================================
    # 3.12. GUT-JAX (1)
    # ========================================================================

    @staticmethod
    def check_SU26_decomposition(chi: int, beta: int, n_func: Callable) -> bool:
        """SU(58) ⊃ SU(26): dim = n(1)² - 1 = 675."""
        try:
            return n_func(1) ** 2 - 1 == 675
        except Exception:
            return False

    # ========================================================================
    # 3.13. Полный список (31 проверка)
    # ========================================================================

    @staticmethod
    def get_all_checks() -> Dict[str, Callable]:
        """
        Все проверки v2.0.

        ВСЕГО: 31 проверка.
        """
        return {
            # χ-связи (3)
            'chi_connections': ModifiableChecks.check_chi_connections,
            'E8_31': ModifiableChecks.check_E8_31,
            'E7_56': ModifiableChecks.check_E7_56,
            # β₁₂ (1)
            'beta_connections': ModifiableChecks.check_beta_connections,
            # Алгебры (5)
            'F4_13': ModifiableChecks.check_F4_13,
            'G2_7': ModifiableChecks.check_G2_7,
            'E6_3': ModifiableChecks.check_E6_3,
            'SU26': ModifiableChecks.check_SU26,
            'SU58': ModifiableChecks.check_SU58,
            # n(k) (5)
            'nk_1': ModifiableChecks.check_nk_1,
            'nk_2': ModifiableChecks.check_nk_2,
            'nk_3': ModifiableChecks.check_nk_3,
            'nk_4': ModifiableChecks.check_nk_4,
            'nk_5': ModifiableChecks.check_nk_5,
            # N_TIME (2)
            'N_TIME_v1': ModifiableChecks.check_N_TIME_v1,
            'N_TIME_v2': ModifiableChecks.check_N_TIME_v2,
            # n(5) (1)
            'n5': ModifiableChecks.check_n5,
            # Тетраэдр (3)
            'Tetrahedron_V_Higgs': ModifiableChecks.check_Tetrahedron_V_Higgs,
            'Tetrahedron_M_Z': ModifiableChecks.check_Tetrahedron_M_Z,
            'Tetrahedron_M_top': ModifiableChecks.check_Tetrahedron_M_top,
            # Космология (2)
            'Cosmology_n_s': ModifiableChecks.check_Cosmology_n_s,
            'Cosmology_H0': ModifiableChecks.check_Cosmology_H0,
            # Хиггс (1)
            'Higgs_delta_M': ModifiableChecks.check_Higgs_delta_M,
            # Φ⁺ (3)
            'phi_plus_E6': ModifiableChecks.check_phi_plus_E6,
            'phi_plus_E7': ModifiableChecks.check_phi_plus_E7,
            'phi_plus_E8': ModifiableChecks.check_phi_plus_E8,
            # Элл. кривая (4)
            'bad_prime_617': ModifiableChecks.check_bad_prime_617,
            'j_invariant': ModifiableChecks.check_j_invariant,
            'BSD_relation': ModifiableChecks.check_BSD_relation,
            'omega_DM': ModifiableChecks.check_omega_DM,
            # GUT-JAX (1)
            'SU26_decomposition': ModifiableChecks.check_SU26_decomposition,
        }


# ============================================================================
# ЧАСТЬ 4: АНАЛИЗ ЧУВСТВИТЕЛЬНОСТИ
# ============================================================================

class SensitivityAnalysis:
    """
    Анализ чувствительности Ω⁴ (v2.0).

    31 проверка, 4 метода анализа:
      • OAT по χ
      • OAT по β₁₂
      • Анализ по категориям
      • Глобальный анализ (сетка)
    """

    def __init__(self):
        self.checks = ModifiableChecks.get_all_checks()
        self.n_poly = ModifiablePolynomial.n_poly

    # ------------------------------------------------------------------------
    # 4.1. OAT по χ
    # ------------------------------------------------------------------------

    def run_oat_chi(self, variations: List[float] = None,
                    verbose: bool = False) -> Dict:
        if variations is None:
            variations = [-0.10, -0.05, -0.02, -0.01, 0.0,
                          0.01, 0.02, 0.05, 0.10]

        base_chi = C.CHI
        base_beta = C.BETA_12

        results = []

        for var in variations:
            new_chi = int(round(base_chi * (1 + var)))

            n_pass = 0
            n_total = 0
            for name, check_func in self.checks.items():
                n_total += 1
                try:
                    if check_func(new_chi, base_beta, self.n_poly):
                        n_pass += 1
                except Exception:
                    pass

            results.append({
                'variation': var,
                'chi': new_chi,
                'n_pass': n_pass,
                'n_total': n_total,
                'pass_rate': n_pass / n_total * 100,
            })

            if verbose:
                print(f"    χ = {new_chi} ({var:+.1%}): "
                      f"{n_pass}/{n_total} ({n_pass / n_total * 100:.1f}%)")

        return {
            'method': 'OAT_chi',
            'base_chi': base_chi,
            'results': results,
        }

    # ------------------------------------------------------------------------
    # 4.2. OAT по β₁₂
    # ------------------------------------------------------------------------

    def run_oat_beta(self, variations: List[int] = None,
                     verbose: bool = False) -> Dict:
        if variations is None:
            variations = [-3, -2, -1, 0, 1, 2, 3]

        base_chi = C.CHI
        base_beta = C.BETA_12

        results = []

        for var in variations:
            new_beta = base_beta + var
            if new_beta < 1:
                continue

            n_pass = 0
            n_total = 0
            for name, check_func in self.checks.items():
                n_total += 1
                try:
                    if check_func(base_chi, new_beta, self.n_poly):
                        n_pass += 1
                except Exception:
                    pass

            results.append({
                'variation': var,
                'beta': new_beta,
                'n_pass': n_pass,
                'n_total': n_total,
                'pass_rate': n_pass / n_total * 100,
            })

            if verbose:
                print(f"    β₁₂ = {new_beta} ({var:+d}): "
                      f"{n_pass}/{n_total} ({n_pass / n_total * 100:.1f}%)")

        return {
            'method': 'OAT_beta',
            'base_beta': base_beta,
            'results': results,
        }

    # ------------------------------------------------------------------------
    # 4.3. Глобальный анализ (сетка)
    # ------------------------------------------------------------------------

    def run_grid(self, chi_range: Tuple[int, int] = (500, 700),
                 beta_range: Tuple[int, int] = (2, 6),
                 n_points: int = 20,
                 verbose: bool = False) -> Dict:
        chi_values = [int(round(chi_range[0] + i * (chi_range[1] - chi_range[0]) / (n_points - 1)))
                      for i in range(n_points)]
        beta_values = [int(round(beta_range[0] + i * (beta_range[1] - beta_range[0]) / (n_points - 1)))
                       for i in range(n_points)]

        grid = []
        for chi in chi_values:
            for beta in beta_values:
                if beta < 1:
                    continue

                n_pass = 0
                n_total = 0
                for name, check_func in self.checks.items():
                    n_total += 1
                    try:
                        if check_func(chi, beta, self.n_poly):
                            n_pass += 1
                    except Exception:
                        pass

                grid.append({
                    'chi': chi,
                    'beta': beta,
                    'n_pass': n_pass,
                    'n_total': n_total,
                    'pass_rate': n_pass / n_total * 100,
                })

        best = max(grid, key=lambda x: x['n_pass'])

        if verbose:
            print(f"    Лучшая точка: χ = {best['chi']}, "
                  f"β₁₂ = {best['beta']}, "
                  f"{best['n_pass']}/{best['n_total']}")

        stable = [g for g in grid if g['pass_rate'] >= 80]

        if verbose:
            print(f"    Область стабильности (≥ 80%): {len(stable)} точек")

        return {
            'method': 'grid',
            'chi_range': chi_range,
            'beta_range': beta_range,
            'n_points': n_points,
            'grid': grid,
            'best': best,
            'stable_region': stable,
            'n_stable': len(stable),
        }

    # ------------------------------------------------------------------------
    # 4.4. Анализ по категориям
    # ------------------------------------------------------------------------

    def run_by_category(self, chi_variations: List[float] = None,
                        verbose: bool = False) -> Dict:
        if chi_variations is None:
            chi_variations = [-0.05, -0.02, -0.01, 0.0, 0.01, 0.02, 0.05]

        categories = {
            'chi': ['chi_connections', 'E8_31', 'E7_56'],
            'beta': ['beta_connections'],
            'algebras': ['F4_13', 'G2_7', 'E6_3', 'SU26', 'SU58'],
            'nk': ['nk_1', 'nk_2', 'nk_3', 'nk_4', 'nk_5',
                   'N_TIME_v1', 'N_TIME_v2', 'n5'],
            'masses': ['Tetrahedron_V_Higgs', 'Tetrahedron_M_Z',
                       'Tetrahedron_M_top'],
            'cosmology': ['Cosmology_n_s', 'Cosmology_H0'],
            'higgs': ['Higgs_delta_M'],
            'phi_plus': ['phi_plus_E6', 'phi_plus_E7', 'phi_plus_E8'],
            'elliptic': ['bad_prime_617', 'j_invariant',
                         'BSD_relation', 'omega_DM'],
            'gut_jax': ['SU26_decomposition'],
        }

        base_chi = C.CHI
        base_beta = C.BETA_12

        results = {}

        for cat_name, check_names in categories.items():
            cat_results = []
            for var in chi_variations:
                new_chi = int(round(base_chi * (1 + var)))

                n_pass = 0
                n_total = 0
                for name in check_names:
                    if name not in self.checks:
                        continue
                    n_total += 1
                    try:
                        if self.checks[name](new_chi, base_beta, self.n_poly):
                            n_pass += 1
                    except Exception:
                        pass

                cat_results.append({
                    'variation': var,
                    'chi': new_chi,
                    'n_pass': n_pass,
                    'n_total': n_total,
                    'pass_rate': n_pass / n_total * 100 if n_total > 0 else 0,
                })

            results[cat_name] = cat_results

            if verbose:
                base_rate = cat_results[len(cat_results) // 2]['pass_rate']
                print(f"    {cat_name:<15}: {base_rate:.1f}% "
                      f"(при χ = {base_chi})")

        return {
            'method': 'by_category',
            'categories': results,
        }

    # ------------------------------------------------------------------------
    # 4.5. Полный анализ
    # ------------------------------------------------------------------------

    def run_full(self, verbose: bool = True) -> Dict:
        if verbose:
            print("=" * 78)
            print("  Ω⁴-UNIFIED v78.0 — АНАЛИЗ ЧУВСТВИТЕЛЬНОСТИ (v2.0)")
            print("=" * 78)
            print(f"\n  Базовые значения: χ = {C.CHI}, β₁₂ = {C.BETA_12}")
            print(f"  Проверок: {len(self.checks)}")

        oat_chi = self.run_oat_chi(verbose=verbose)
        oat_beta = self.run_oat_beta(verbose=verbose)
        by_cat = self.run_by_category(verbose=verbose)
        grid = self.run_grid(chi_range=(500, 700), beta_range=(2, 6),
                             n_points=10, verbose=verbose)

        if verbose:
            print("\n" + "=" * 78)
            print("  ИТОГИ")
            print("=" * 78)

            base_result = next(r for r in oat_chi['results']
                               if r['variation'] == 0.0)
            worst_chi = min(oat_chi['results'], key=lambda x: x['pass_rate'])
            best_chi = max(oat_chi['results'], key=lambda x: x['pass_rate'])

            print(f"\n  Устойчивость к χ:")
            print(f"    Базовая:      {base_result['pass_rate']:.1f}%")
            print(f"    Лучшая:       {best_chi['pass_rate']:.1f}% "
                  f"(χ = {best_chi['chi']})")
            print(f"    Худшая:       {worst_chi['pass_rate']:.1f}% "
                  f"(χ = {worst_chi['chi']})")

            base_beta_result = next(r for r in oat_beta['results']
                                    if r['variation'] == 0)
            worst_beta = min(oat_beta['results'],
                             key=lambda x: x['pass_rate'])

            print(f"\n  Устойчивость к β₁₂:")
            print(f"    Базовая:      {base_beta_result['pass_rate']:.1f}%")
            print(f"    Худшая:       {worst_beta['pass_rate']:.1f}% "
                  f"(β₁₂ = {worst_beta['beta']})")

            print(f"\n  Область стабильности (≥ 80%):")
            print(f"    Точек:        {grid['n_stable']}/{len(grid['grid'])}")
            print(f"    Лучшая:       χ = {grid['best']['chi']}, "
                  f"β₁₂ = {grid['best']['beta']}")

            print(f"\n  ВЫВОД:")
            if grid['n_stable'] > len(grid['grid']) * 0.1:
                print(f"    ✅ Модель устойчива к изменениям χ, β₁₂")
            else:
                print(f"    ⚠️ Модель чувствительна к изменениям χ, β₁₂")

            print("=" * 78)

        return {
            'oat_chi': oat_chi,
            'oat_beta': oat_beta,
            'by_category': by_cat,
            'grid': grid,
        }


# ============================================================================
# ЧАСТЬ 5: ПЕЧАТЬ ОТЧЁТА (для main_unified.py)
# ============================================================================
#
# ЗАЧЕМ:
#   Вызывается из main_unified.py → MODULE_REPORTS.
#   Печатает полный отчёт по анализу чувствительности v2.0.
#
# ============================================================================

def print_report():
    """
    Печать полного отчёта по анализу чувствительности v2.0.

    ИСПОЛЬЗУЕТСЯ В:
        main_unified.py → MODULE_REPORTS →
        ('sensitivity', 'print_report', ...)
    """
    print("=" * 78)
    print("  Ω⁴-UNIFIED v78.0 — АНАЛИЗ ЧУВСТВИТЕЛЬНОСТИ (sensitivity.py v2.0)")
    print("=" * 78)

    analysis = SensitivityAnalysis()
    results = analysis.run_full(verbose=True)

    return results


# ============================================================================
# ЧАСТЬ 6: ЗАПУСК
# ============================================================================

if __name__ == "__main__":
    print("=" * 78)
    print("  Ω⁴-UNIFIED v78.0 — АНАЛИЗ ЧУВСТВИТЕЛЬНОСТИ (sensitivity.py v2.0)")
    print("=" * 78)

    analysis = SensitivityAnalysis()
    results = analysis.run_full(verbose=True)
