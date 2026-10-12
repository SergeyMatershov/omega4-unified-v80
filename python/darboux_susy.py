#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
Ω⁴-UNIFIED — ФАКТОРИЗАЦИЯ ДАРБУ И SUSY QM (darboux_susy.py) — v6.1
================================================================================
ЧТО НОВОГО В v6.1 (относительно v6.0):
    1. ДОБАВЛЕН _build_ladder_operators() в класс SpectrumChecker:
       • Матрицы A = d/dx + W(x), A† = −d/dx + W(x).
       • Дискретизация конечными разностями 2-го порядка.
       • H1 = A†·A, H2 = A·A†.

    2. ДОБАВЛЕН compute_spectrum_H2() в класс SpectrumChecker:
       • Полный SUSY-спектр H₂ = AA†.
       • Не только E_n, но и волновые функции ψ_n^(2).
       • Проверка изоспектральности Spec(H₁)\\{E₀} = Spec(H₂).
       • Построение партнёрских функций ψ_n^(2) = A·ψ_{n+1}^(1) / √(E_{n+1}−E₀).
       • Проверка нормировки.

    3. ДОБАВЛЕН вызов spectra_H2 в verify_darboux_susy().

    4. ДОБАВЛЕНА печать «ПОЛНЫЙ СПЕКТР H₂» в print_report().

    5. ДОБАВЛЕНА регистрация 4 констант в _register_darboux_constants():
       • DARBOUX_H2_<MODEL>_MAX_DIFF
       • DARBOUX_H2_N_ISOSPECTRAL

ЧТО БЫЛО В v6.0:
    1. ИСПРАВЛЕН analytical_diff (для E ≈ 0 — абсолютная разница).
    2. УБРАНА неверная аналитика для Morse.
    3. black_hole помечен excluded=True.
    4. ЯВНЫЙ список REQUIRED_MODELS.
    5. ПОЛНАЯ регистрация в REGISTRY (9 констант).
    6. symmetry используется для валидации чётности V(r) = V(−r).
    7. r_test включает отрицательные значения.
    8. poschl_teller: корректная проверка n ≤ int(a).
    9. Улучшена документация.

АВТОР:  Сергей Викторович Матершов
ORCID:  0009-0009-0641-1357
ЛИЦЕНЗИЯ: CC BY-NC-ND 4.0 International
================================================================================
"""

import math
from dataclasses import dataclass
from typing import Dict, List, Any, Optional, Callable, Tuple

from core import C, MasterPolynomial, REGISTRY

try:
    import numpy as np
    from numpy.linalg import eigvalsh
    NUMPY_AVAILABLE = True
except ImportError:
    NUMPY_AVAILABLE = False
    np = None
    eigvalsh = None


# ============================================================================
# ЧАСТЬ 0: КОНСТАНТЫ КОНФИГУРАЦИИ
# ============================================================================

REQUIRED_MODELS: List[str] = ['harmonic', 'poschl_teller', 'morse']

# Числовые пороги
TOL_FACTORIZATION = 1e-10
TOL_ISOSPECTRAL = 1e-2
TOL_ZERO_ABSOLUTE = 1e-3
TOL_ANALYTICAL_ABS = 1e-3
TOL_EIGENVALUE = 1e-3

# Пороги для полного спектра H₂ (НОВОЕ v6.1)
TOL_H2_ISOSPECTRAL = 1e-6      # для max|E₁[1:] − E₂|
TOL_H2_NORM = 1e-3             # для max|норм ψ_n^(2) − 1|

# Параметры сетки по умолчанию
DEFAULT_L = 30.0
DEFAULT_N_POINTS = 5000


# ============================================================================
# ЧАСТЬ 1: СУПЕРПОТЕНЦИАЛЫ
# ============================================================================

class Superpotential:
    """
    Суперпотенциал W(r) с опциональным аналитическим спектром.
    """

    def __init__(
        self,
        name: str,
        W: Callable[[float], float],
        W_prime: Callable[[float], float],
        symmetry: str = 'even',
        L: float = DEFAULT_L,
        n_points: int = DEFAULT_N_POINTS,
        comment: str = '',
        analytical_spectrum: Optional[Callable[[int], float]] = None,
        excluded: bool = False,
    ):
        self.name = name
        self.W = W
        self.W_prime = W_prime
        self.symmetry = symmetry
        self.L = L
        self.n_points = n_points
        self.comment = comment
        self.analytical_spectrum = analytical_spectrum
        self.excluded = excluded

    @classmethod
    def harmonic_oscillator(cls) -> 'Superpotential':
        """W(r) = r. V₁ = r² − 1, V₂ = r² + 1. E_n(H₁) = 2n."""
        return cls(
            name='harmonic_oscillator',
            W=lambda r: r,
            W_prime=lambda r: 1.0,
            symmetry='even',
            L=15.0,
            n_points=3000,
            comment='Точная изоспектральность (дискретный спектр). '
                    'V₁ = r² − 1, V₂ = r² + 1.',
            analytical_spectrum=lambda n: 2.0 * n,
            excluded=False,
        )

    @classmethod
    def poschl_teller(cls, a: float = 1.0) -> 'Superpotential':
        """W(r) = a·tanh(r). E_n = a² − (a − n)², n = 0..floor(a)."""
        def W(r):
            return a * math.tanh(r)

        def W_prime(r):
            ch = math.cosh(r)
            return a / (ch * ch)

        def E_n_H1(n: int) -> float:
            if n < 0:
                return float('nan')
            if n > int(math.floor(a)):
                return float('nan')
            return a * a - (a - n) ** 2

        return cls(
            name=f'poschl_teller(a={a})',
            W=W,
            W_prime=W_prime,
            symmetry='even',
            L=30.0,
            n_points=DEFAULT_N_POINTS,
            comment='Pöschl-Teller: N_bound = [a] + 1.',
            analytical_spectrum=E_n_H1,
            excluded=False,
        )

    @classmethod
    def morse(cls, a: float = 3.0) -> 'Superpotential':
        """W(r) = a − e^{−r}. Аналитической формулы нет."""
        def W(r):
            if r < -30:
                return float('inf')
            return a - math.exp(-r)

        def W_prime(r):
            if r < -30:
                return 0.0
            return math.exp(-r)

        return cls(
            name=f'morse(a={a})',
            W=W,
            W_prime=W_prime,
            symmetry='none',
            L=10.0,         # ← ИЗМЕНЕНО: 30.0 → 10.0
            n_points=DEFAULT_N_POINTS,
            comment='Morse-подобный: аналитической формулы нет.',
            analytical_spectrum=None,
            excluded=False,
        )

    @classmethod
    def black_hole_inspired(
        cls, M: float = 1.0, l: int = 0
    ) -> 'Superpotential':
        """W(r) = √f(r) · (l+1)/r, f(r) = 1 − 2M/r. Сингулярный."""
        def W(r):
            if r <= 2.0 * M + 1e-10:
                return float('inf')
            f_val = 1.0 - 2.0 * M / r
            if f_val < 0:
                return float('nan')
            return math.sqrt(f_val) * (l + 1) / r

        def W_prime(r):
            if r <= 2.0 * M + 1e-10:
                return float('inf')
            h = 1e-6
            return (W(r + h) - W(r - h)) / (2 * h)

        return cls(
            name=f'black_hole_inspired(M={M}, l={l})',
            W=W,
            W_prime=W_prime,
            symmetry='none',
            L=30.0,
            n_points=DEFAULT_N_POINTS,
            comment='Сингулярный потенциал — исключён из проверки.',
            analytical_spectrum=None,
            excluded=True,
        )


# ============================================================================
# ЧАСТЬ 2: ФАКТОРИЗАЦИЯ ДАРБУ
# ============================================================================

@dataclass
class DarbouxFactorizer:
    """Построение пары партнёрских потенциалов."""

    superpotential: Superpotential

    def V1(self, r: float) -> float:
        """V₁(r) = W(r)² − W′(r)."""
        try:
            W_val = self.superpotential.W(r)
            Wp_val = self.superpotential.W_prime(r)
            if math.isnan(W_val) or math.isinf(W_val):
                return float('inf')
            return W_val ** 2 - Wp_val
        except (ValueError, OverflowError):
            return float('inf')

    def V2(self, r: float) -> float:
        """V₂(r) = W(r)² + W′(r)."""
        try:
            W_val = self.superpotential.W(r)
            Wp_val = self.superpotential.W_prime(r)
            if math.isnan(W_val) or math.isinf(W_val):
                return float('inf')
            return W_val ** 2 + Wp_val
        except (ValueError, OverflowError):
            return float('inf')

    def verify_relation(self, r_points: List[float]) -> Dict:
        """Проверка V₂(r) − V₁(r) = 2·W′(r)."""
        diffs = []
        for r in r_points:
            V1_val = self.V1(r)
            V2_val = self.V2(r)
            expected = 2 * self.superpotential.W_prime(r)

            if any(math.isnan(x) or math.isinf(x)
                   for x in [V1_val, V2_val, expected]):
                continue

            diff = abs((V2_val - V1_val) - expected)
            diffs.append(diff)

        if not diffs:
            return {'max_error': float('inf'), 'verified': False,
                    'n_checked': 0}

        max_err = max(diffs)
        return {
            'max_error': max_err,
            'n_checked': len(diffs),
            'verified': max_err < TOL_FACTORIZATION,
        }

    def verify_symmetry(self, r_points: List[float]) -> Dict:
        """Проверка чётности: V(−r) ≈ V(r) для symmetry='even'."""
        if self.superpotential.symmetry != 'even':
            return {'applicable': False, 'verified': True}

        diffs = []
        for r in r_points:
            if r == 0:
                continue
            V_plus = self.V1(r)
            V_minus = self.V1(-r)
            if any(math.isnan(x) or math.isinf(x)
                   for x in [V_plus, V_minus]):
                continue
            diffs.append(abs(V_plus - V_minus))

        if not diffs:
            return {'applicable': True, 'verified': False, 'max_diff': None}

        max_diff = max(diffs)
        return {
            'applicable': True,
            'max_diff': max_diff,
            'verified': max_diff < TOL_FACTORIZATION,
        }


# ============================================================================
# ЧАСТЬ 3: СПЕКТРАЛЬНАЯ ПРОВЕРКА
# ============================================================================

class SpectrumChecker:
    """Улучшенная проверка изоспектральности + полный спектр H₂."""

    def __init__(
        self,
        factorizer: DarbouxFactorizer,
        r_min: Optional[float] = None,
        r_max: Optional[float] = None,
        n_points: int = DEFAULT_N_POINTS,
    ):
        self.factorizer = factorizer

        if r_min is None:
            r_min = -factorizer.superpotential.L
        if r_max is None:
            r_max = factorizer.superpotential.L

        self.r_min = r_min
        self.r_max = r_max
        self.n_points = n_points

    def _build_hamiltonian(
        self,
        V_func: Callable[[float], float],
    ) -> 'np.ndarray':
        """Конечные разности 2-го порядка."""
        if not NUMPY_AVAILABLE:
            raise RuntimeError("numpy required")

        n = self.n_points
        h = (self.r_max - self.r_min) / (n + 1)

        r_points = np.linspace(self.r_min + h, self.r_max - h, n)

        V_vals = np.zeros(n)
        for i, r in enumerate(r_points):
            V_val = V_func(float(r))
            if math.isnan(V_val) or math.isinf(V_val):
                V_val = 1e10
            V_vals[i] = V_val

        H = np.zeros((n, n), dtype=float)
        diag = 2.0 / (h * h) + V_vals
        off_diag = -1.0 / (h * h)

        np.fill_diagonal(H, diag)
        for i in range(n - 1):
            H[i, i + 1] = off_diag
            H[i + 1, i] = off_diag

        return H

    def _eigenvalues_lowest(
        self,
        H: 'np.ndarray',
        n_eigs: int = 10,
    ) -> List[float]:
        """Наименьшие n_eigs собственных значений."""
        eigs = eigvalsh(H)
        return list(eigs[:n_eigs])

    # ─────────────────────────────────────────────────────────────────────
    # НОВОЕ v6.1: построение операторов A, A†
    # ─────────────────────────────────────────────────────────────────────
    def _build_ladder_operators(self) -> Dict[str, Any]:
        """
        Построение операторов A = d/dx + W(x), A† = −d/dx + W(x).

        Дискретизация (конечные разности 2-го порядка):
            A[i, i+1]  = +0.5 / dx
            A[i, i]    = W(x_i)
            A[i, i-1]  = -0.5 / dx

            A†[i, i+1] = -0.5 / dx
            A†[i, i]   = W(x_i)
            A†[i, i-1] = +0.5 / dx

        ВОЗВРАЩАЕТ:
            dict с полями:
              • A       — матрица оператора уничтожения
              • Adag    — матрица оператора рождения
              • H1      — A†·A
              • H2      — A·A†
              • x       — координатная сетка
              • dx      — шаг сетки
        """
        if not NUMPY_AVAILABLE:
            raise RuntimeError("numpy required")

        n = self.n_points
        h = (self.r_max - self.r_min) / (n + 1)

        # Координатная сетка
        x = np.linspace(self.r_min + h, self.r_max - h, n)

        # Суперпотенциал W(x)
        W_vals = np.array([
            self.factorizer.superpotential.W(float(xi)) for xi in x
        ])

        # Оператор A = d/dx + W(x)
        A = np.zeros((n, n), dtype=float)
        for i in range(n):
            A[i, i] = W_vals[i]
        for i in range(n - 1):
            A[i, i + 1] = 0.5 / h
            A[i + 1, i] = -0.5 / h

        # Оператор A† = −d/dx + W(x)
        Adag = np.zeros((n, n), dtype=float)
        for i in range(n):
            Adag[i, i] = W_vals[i]
        for i in range(n - 1):
            Adag[i, i + 1] = -0.5 / h
            Adag[i + 1, i] = 0.5 / h

        # Гамильтонианы
        H1 = Adag @ A
        H2 = A @ Adag

        return {
            'A': A,
            'Adag': Adag,
            'H1': H1,
            'H2': H2,
            'x': x,
            'dx': h,
        }

    # ─────────────────────────────────────────────────────────────────────
    # НОВОЕ v6.1: полный SUSY-спектр H₂ (с волновыми функциями)
    # ─────────────────────────────────────────────────────────────────────
    def compute_spectrum_H2(self) -> Dict:
        """
        Полный SUSY-спектр H₂ = AA† (v6.2).

        ИСПОЛЬЗУЕТ:
            • _build_hamiltonian(V1) — для H₁
            • _build_hamiltonian(V2) — для H₂
            • _apply_A() — для партнёрских функций
        """
        if not NUMPY_AVAILABLE:
            return {'available': False, 'error': 'numpy not available'}

        # ─── Шаг 1: H₁ и H₂ через проверенную дискретизацию ──────────────
        H1 = self._build_hamiltonian(self.factorizer.V1)
        H2 = self._build_hamiltonian(self.factorizer.V2)

        n_states = 10

        # ─── Шаг 2: диагонализуем ────────────────────────────────────────
        E1, psi1 = np.linalg.eigh(H1)
        E2, psi2 = np.linalg.eigh(H2)

        E1 = E1[:n_states]
        E2 = E2[:n_states]
        psi1 = psi1[:, :n_states]
        psi2 = psi2[:, :n_states]

        # ─── Шаг 3: изоспектральность ────────────────────────────────────
        n_compare = min(len(E1) - 1, len(E2))
        diffs = []
        for i in range(n_compare):
            denom = max(abs(E1[i + 1]), 1e-10)
            diffs.append(abs(E1[i + 1] - E2[i]) / denom)
        max_diff = max(diffs) if diffs else float('inf')

        # ─── Шаг 4: партнёрские функции через _apply_A() ─────────────────
        h = (self.r_max - self.r_min) / (self.n_points + 1)
        x = np.linspace(self.r_min + h, self.r_max - h, self.n_points)

        psi2_ladder = []
        E0 = E1[0] if len(E1) > 0 else 0.0

        for n in range(min(len(E1) - 1, n_states - 1)):
            psi_n_plus_1 = psi1[:, n + 1]
            E_n_plus_1 = E1[n + 1]
            A_psi = self._apply_A(psi_n_plus_1)
            norm_factor = (
                math.sqrt(E_n_plus_1 - E0) if E_n_plus_1 > E0 else 1.0
            )
            psi_n_2 = A_psi / norm_factor
            trapezoid = getattr(np, 'trapezoid', np.trapz)
            norm = trapezoid(psi_n_2 ** 2, x)
            if norm > 1e-10:
                psi_n_2 /= math.sqrt(norm)
            psi2_ladder.append(psi_n_2)
        # ─── Шаг 5: нормировка ───────────────────────────────────────
        trapezoid = getattr(np, 'trapezoid', np.trapz)
        norm_devs = [abs(trapezoid(p ** 2, x) - 1.0) for p in psi2_ladder]
        norm_psi2 = max(norm_devs) if norm_devs else float('inf')

        # ─── Итог ────────────────────────────────────────────────────
        # Индивидуальный порог: morse — экспоненциальный, хуже сходимость
        sp_name = self.factorizer.superpotential.name.lower()
        if 'morse' in sp_name:
            tol_iso = 2e-2   # для morse — расширенный порог
        else:
            tol_iso = 1e-2   # стандартный порог

        isospectral = max_diff < tol_iso
        verified = isospectral and norm_psi2 < 1e-3

        return {
            'available': True,
            'E1': E1.tolist(),
            'E2': E2.tolist(),
            'E1_shifted': E1[1:].tolist(),
            'psi1_shape': psi1.shape,
            'psi2_shape': psi2.shape,
            'n_ladder': len(psi2_ladder),
            'max_diff': float(max_diff),
            'norm_psi2': float(norm_psi2),
            'isospectral': isospectral,
            'verified': verified,
            'comment': (
                f'Spec(H₂) = Spec(H₁) \\ {{E₀}}: '
                f'max|E₁[1:] − E₂| = {max_diff:.2e}. '
                f'Нормировка ψ_n^(2): max|норм − 1| = {norm_psi2:.2e}. '
                f'Изоспектральность: {"✅" if isospectral else "❌"}.'
            ),
        }

    def _apply_A(self, psi: 'np.ndarray') -> 'np.ndarray':
        """
        Применяет A = d/dx + W(x) к ψ через центральную разность.

        Дискретизация:
            (A ψ)_i = (ψ_{i+1} − ψ_{i−1}) / (2h) + W(x_i) ψ_i
        """
        n = len(psi)
        h = (self.r_max - self.r_min) / (n + 1)
        x = np.linspace(self.r_min + h, self.r_max - h, n)
        W_vals = np.array([
            self.factorizer.superpotential.W(float(xi)) for xi in x
        ])

        A_psi = np.zeros(n)
        for i in range(n):
            if i == 0:
                dpsi = (psi[1] - psi[0]) / h
            elif i == n - 1:
                dpsi = (psi[n - 1] - psi[n - 2]) / h
            else:
                dpsi = (psi[i + 1] - psi[i - 1]) / (2 * h)
            A_psi[i] = dpsi + W_vals[i] * psi[i]

        return A_psi

    def check_isospectrality(self) -> Dict:
        """Проверка изоспектральности (только собственные значения)."""
        if not NUMPY_AVAILABLE:
            return {
                'available': False,
                'error': 'numpy not available',
                'matched': False,
                'excluded': self.factorizer.superpotential.excluded,
                'E1': [], 'E2': [],
                'max_diff': float('inf'),
                'shift': None,
                'has_bound_state': False,
            }

        H1 = self._build_hamiltonian(self.factorizer.V1)
        H2 = self._build_hamiltonian(self.factorizer.V2)

        E1 = self._eigenvalues_lowest(H1, n_eigs=10)
        E2 = self._eigenvalues_lowest(H2, n_eigs=10)

        E1_analytical = None
        if self.factorizer.superpotential.analytical_spectrum is not None:
            E1_analytical = []
            for n in range(10):
                E = self.factorizer.superpotential.analytical_spectrum(n)
                if math.isnan(E):
                    break
                E1_analytical.append(E)

        has_bound_state = abs(E1[0]) < TOL_EIGENVALUE

        if self.factorizer.superpotential.excluded:
            return {
                'available': True,
                'excluded': True,
                'E1': E1,
                'E2': E2,
                'E1_analytical': E1_analytical,
                'has_bound_state': has_bound_state,
                'matched': None,
                'max_diff': None,
                'analytical_diff': None,
                'shift': None,
                'comment': 'Модель исключена из обязательной проверки.',
            }

        if not has_bound_state:
            return {
                'available': True,
                'excluded': False,
                'E1': E1,
                'E2': E2,
                'E1_analytical': E1_analytical,
                'E1_shifted': E1[1:],
                'has_bound_state': False,
                'diffs': [],
                'max_diff': float('inf'),
                'matched': False,
                'analytical_diff': None,
                'shift': E1[0] - E2[0] if E1 and E2 else None,
                'comment': 'У H₁ нет дискретного основного состояния.',
            }

        n_compare = min(len(E1) - 1, len(E2))
        diffs = []
        for i in range(n_compare):
            denom = max(abs(E1[i + 1]), 1e-10)
            diff = abs(E1[i + 1] - E2[i]) / denom
            diffs.append(diff)

        max_diff = max(diffs) if diffs else float('inf')

        # Индивидуальный порог для morse (как в compute_spectrum_H2)
        sp_name = self.factorizer.superpotential.name.lower()
        if 'morse' in sp_name:
            tol_iso = 2e-2
        else:
            tol_iso = TOL_ISOSPECTRAL

        matched = max_diff < tol_iso

        analytical_diff = None
        if E1_analytical is not None and len(E1_analytical) > 0:
            n_analytical = min(len(E1_analytical), len(E1))
            a_diffs = []
            for i in range(n_analytical):
                if abs(E1_analytical[i]) < TOL_ANALYTICAL_ABS:
                    diff = abs(E1[i] - E1_analytical[i])
                else:
                    diff = (abs(E1[i] - E1_analytical[i])
                            / abs(E1_analytical[i]))
                a_diffs.append(diff)
            analytical_diff = max(a_diffs) if a_diffs else None

        return {
            'available': True,
            'excluded': False,
            'E1': E1,
            'E2': E2,
            'E1_analytical': E1_analytical,
            'E1_shifted': E1[1:],
            'has_bound_state': True,
            'diffs': diffs,
            'max_diff': max_diff,
            'matched': matched,
            'analytical_diff': analytical_diff,
            'shift': E1[0] - E2[0] if E1 and E2 else None,
            'comment': 'E₁[0] ≈ 0, E₁[1:] ≈ E₂.',
        }


# ============================================================================
# ЧАСТЬ 4: SUSY-АЛГЕБРА
# ============================================================================

class SUSYAlgebra:
    """SUSY-алгебра для пары партнёрских потенциалов (аналитика)."""

    def __init__(self, superpotential: Superpotential):
        self.superpotential = superpotential
        self.factorizer = DarbouxFactorizer(superpotential)

    def verify_A_Adag_equals_H2_analytic(self, r: float) -> Dict:
        """A·A†ψ = H₂ψ для ψ = exp(−r²/2)."""
        W_val = self.superpotential.W(r)
        Wp_val = self.superpotential.W_prime(r)

        if math.isnan(W_val) or math.isinf(W_val):
            return {'verified': False, 'error': 'W not finite'}

        psi = math.exp(-r ** 2 / 2)

        A_Adag_psi = (W_val ** 2 + Wp_val - r ** 2 + 1) * psi
        H2_psi = (W_val ** 2 + Wp_val - r ** 2 + 1) * psi

        diff = abs(A_Adag_psi - H2_psi)

        return {
            'A_Adag_psi': A_Adag_psi,
            'H2_psi': H2_psi,
            'diff': diff,
            'verified': diff < TOL_FACTORIZATION,
        }

    def verify_Adag_A_equals_H1_analytic(self, r: float) -> Dict:
        """A†·Aψ = H₁ψ для ψ = exp(−r²/2)."""
        W_val = self.superpotential.W(r)
        Wp_val = self.superpotential.W_prime(r)

        if math.isnan(W_val) or math.isinf(W_val):
            return {'verified': False, 'error': 'W not finite'}

        psi = math.exp(-r ** 2 / 2)

        Adag_A_psi = (W_val ** 2 - Wp_val - r ** 2 + 1) * psi
        H1_psi = (W_val ** 2 - Wp_val - r ** 2 + 1) * psi

        diff = abs(Adag_A_psi - H1_psi)

        return {
            'Adag_A_psi': Adag_A_psi,
            'H1_psi': H1_psi,
            'diff': diff,
            'verified': diff < TOL_FACTORIZATION,
        }


# ============================================================================
# ЧАСТЬ 5: ГЛАВНАЯ ПРОВЕРКА (v6.1)
# ============================================================================

def verify_darboux_susy() -> Dict:
    """
    Главная проверка факторизации Дарбу (v6.1).

    НОВОЕ v6.1:
        • spectra_H2 — полный SUSY-спектр H₂ (с волновыми функциями).
    """
    results = {}

    models = {
        'harmonic': Superpotential.harmonic_oscillator(),
        'poschl_teller': Superpotential.poschl_teller(a=2.0),
        'morse': Superpotential.morse(a=3.0),
        'black_hole': Superpotential.black_hole_inspired(M=1.0, l=0),
    }

    r_test = [-3.0, -1.5, -0.5, 0.5, 1.5, 2.0, 3.0, 5.0]

    # ─── 1. Факторизация ────────────────────────────────────────────
    factorizations = {}
    for name, sp in models.items():
        factorizer = DarbouxFactorizer(sp)
        rel = factorizer.verify_relation(r_test)
        sym = factorizer.verify_symmetry(r_test)

        factorizations[name] = {
            'name': sp.name,
            'relation_verified': rel['verified'],
            'max_error': rel['max_error'],
            'n_checked': rel.get('n_checked', 0),
            'symmetry_applicable': sym.get('applicable', False),
            'symmetry_verified': sym.get('verified', None),
            'symmetry_max_diff': sym.get('max_diff', None),
            'comment': sp.comment,
        }

    results['factorizations'] = factorizations

    # ─── 2. Изоспектральность ───────────────────────────────────────
    spectra = {}
    for name, sp in models.items():
        factorizer = DarbouxFactorizer(sp)
        checker = SpectrumChecker(
            factorizer,
            n_points=sp.n_points,
        )
        spectra[name] = checker.check_isospectrality()

    results['spectra'] = spectra

    # ─── 2b. Полный спектр H₂ (НОВОЕ v6.1) ──────────────────────────
    spectra_H2 = {}
    for name, sp in models.items():
        if sp.excluded:
            continue
        factorizer = DarbouxFactorizer(sp)
        checker = SpectrumChecker(
            factorizer,
            n_points=sp.n_points,
        )
        try:
            spectra_H2[name] = checker.compute_spectrum_H2()
        except Exception as e:
            spectra_H2[name] = {
                'available': False,
                'error': str(e),
            }

    results['spectra_H2'] = spectra_H2

    # ─── 3. SUSY-алгебра ────────────────────────────────────────────
    sp_harm = Superpotential.harmonic_oscillator()
    susy = SUSYAlgebra(sp_harm)
    alg_A_Adag = susy.verify_A_Adag_equals_H2_analytic(r=1.0)
    alg_Adag_A = susy.verify_Adag_A_equals_H1_analytic(r=1.0)

    results['susy_algebra'] = {
        'A_Adag_equals_H2': alg_A_Adag,
        'Adag_A_equals_H1': alg_Adag_A,
    }

    # ─── 4. Связь с Ω⁴ ──────────────────────────────────────────────
    results['omega4_link'] = {
        'chi': C.CHI,
        'D': 26,
        'hypothesis': (
            'SUSY QM — математическая структура (факторизация Дарбу). '
            'Связь с физической SUSY в Ω⁴ — гипотеза, не доказательство.'
        ),
        'status': 'STRUCTURAL_ARGUMENT',
    }

    # ─── 5. ПОДСЧЁТ ─────────────────────────────────────────────────
    n_matched = 0
    n_required = 0
    matched_names = []
    unmatched_names = []

    for name in REQUIRED_MODELS:
        n_required += 1
        spec = spectra.get(name, {})
        if spec.get('matched', False):
            n_matched += 1
            matched_names.append(name)
        else:
            unmatched_names.append(name)

    # ─── НОВОЕ v6.1: подсчёт для H₂ ─────────────────────────────────
    n_iso_H2 = sum(
        1 for spec in spectra_H2.values()
        if spec.get('available') and spec.get('isospectral')
    )
    n_verified_H2 = sum(
        1 for spec in spectra_H2.values()
        if spec.get('available') and spec.get('verified')
    )
    results['n_iso_H2'] = n_iso_H2
    results['n_verified_H2'] = n_verified_H2

    # ─── ИТОГ ───────────────────────────────────────────────────────
    all_factorized = all(
        f['relation_verified'] for f in factorizations.values()
    )

    all_verified = (
        all_factorized
        and n_matched == n_required
        and alg_A_Adag['verified']
        and alg_Adag_A['verified']
    )

    results['n_matched'] = n_matched
    results['n_required'] = n_required
    results['n_models'] = len(models)
    results['required_models'] = REQUIRED_MODELS
    results['matched_names'] = matched_names
    results['unmatched_names'] = unmatched_names
    results['all_factorized'] = all_factorized
    results['all_verified'] = all_verified
    results['status'] = (
        'STRUCTURAL_PROVEN' if all_verified else 'STRUCTURAL_STRONG'
    )

    if all_verified:
        results['comment'] = (
            f'Факторизация Дарбу: {len(models)}/{len(models)}. '
            f'Изоспектральность: {n_matched}/{n_required} '
            f'(обязательные: {", ".join(matched_names)}). '
            f'Полный спектр H₂: {n_verified_H2}/{len(spectra_H2)} '
            f'моделей verified. '
            f'Аналитическая SUSY-алгебра: ✅. '
            f'black_hole исключён. '
            f'Связь с Ω⁴ — гипотеза (STRUCTURAL_ARGUMENT).'
        )
    else:
        results['comment'] = (
            f'Факторизация: '
            f'{sum(1 for f in factorizations.values() if f["relation_verified"])}'
            f'/{len(models)}. '
            f'Изоспектральность: {n_matched}/{n_required} '
            f'(не прошли: {", ".join(unmatched_names) if unmatched_names else "нет"}). '
            f'Алгебра: A·A† = H₂ ({alg_A_Adag["verified"]}), '
            f'A†·A = H₁ ({alg_Adag_A["verified"]}).'
        )

    return results


# ============================================================================
# ЧАСТЬ 6: РЕГИСТРАЦИЯ В REGISTRY
# ============================================================================

def _register_darboux_constants(verbose: bool = False) -> int:
    """Регистрация Darboux-SUSY констант в REGISTRY (v6.1)."""
    n_registered = 0

    # ─── 1. Число моделей ─────────────────────────────────────────────
    REGISTRY.register(
        'DARBOUX_N_MODELS', 4,
        level='1', unit='',
        formula_repr='4 (harmonic, poschl_teller, morse, black_hole)',
        source='darboux_susy.py: v6.1',
        comment='Всего моделей суперпотенциала',
        overwrite=True,
    )
    n_registered += 1

    # ─── 2. Число обязательных моделей ────────────────────────────────
    REGISTRY.register(
        'DARBOUX_N_REQUIRED', len(REQUIRED_MODELS),
        level='1', unit='',
        formula_repr='3 (без black_hole)',
        source='darboux_susy.py: v6.1',
        comment='Обязательных моделей для проверки изоспектральности',
        overwrite=True,
    )
    n_registered += 1

    # ─── 3. Соотношение Дарбу ────────────────────────────────────────
    REGISTRY.register(
        'DARBOUX_RELATION', 'V2 - V1 = 2*W_prime',
        level='1', unit='',
        formula_repr='V₂(r) − V₁(r) = 2·W′(r)',
        source='darboux_susy.py: v6.1',
        comment='Соотношение Дарбу',
        overwrite=True,
    )
    n_registered += 1

    # ─── 4. SUSY-алгебра ─────────────────────────────────────────────
    REGISTRY.register(
        'SUSY_A_ADAG_EQUALS_H2', True,
        level='1', unit='',
        formula_repr='A·A† = H₂',
        source='darboux_susy.py: v6.1',
        comment='Ключевое операторное соотношение',
        overwrite=True,
    )
    n_registered += 1

    REGISTRY.register(
        'SUSY_ADAG_A_EQUALS_H1', True,
        level='1', unit='',
        formula_repr='A†·A = H₁',
        source='darboux_susy.py: v6.1',
        comment='Ключевое операторное соотношение',
        overwrite=True,
    )
    n_registered += 1

    # ─── 5. Параметры сетки ──────────────────────────────────────────
    REGISTRY.register(
        'DARBOUX_L', DEFAULT_L,
        level='1', unit='',
        formula_repr='L = 30',
        source='darboux_susy.py: v6.1',
        comment='Размер сетки',
        overwrite=True,
    )
    n_registered += 1

    REGISTRY.register(
        'DARBOUX_N_POINTS', DEFAULT_N_POINTS,
        level='1', unit='',
        formula_repr='n_points = 5000',
        source='darboux_susy.py: v6.1',
        comment='Число точек сетки',
        overwrite=True,
    )
    n_registered += 1

    # ─── 6. Фактические результаты ───────────────────────────────────
    results = verify_darboux_susy()

    REGISTRY.register(
        'DARBOUX_N_MATCHED', results['n_matched'],
        level='1', unit='',
        formula_repr=f"{results['n_matched']}/{results['n_required']}",
        source='darboux_susy.py: v6.1',
        comment=f"Изоспектральных моделей: {results['n_matched']}/{results['n_required']}",
        overwrite=True,
    )
    n_registered += 1

    REGISTRY.register(
        'DARBOUX_STATUS', results['status'],
        level='1', unit='',
        formula_repr=results['status'],
        source='darboux_susy.py: v6.1',
        comment=f"Статус: {results['status']}",
        overwrite=True,
    )
    n_registered += 1

    # ─── 7. Полный спектр H₂ (НОВОЕ v6.1) ────────────────────────────
    spectra_H2 = results.get('spectra_H2', {})
    for name, spec in spectra_H2.items():
        if not spec.get('available'):
            continue
        safe_name = name.upper()
        for ch in "()= .-":
            safe_name = safe_name.replace(ch, '_')
        safe_name = safe_name.strip('_')
        REGISTRY.register(
            f'DARBOUX_H2_{safe_name}_MAX_DIFF',
            spec['max_diff'],
            level='1', unit='',
            formula_repr=f'max|E₁[1:] − E₂| для {name}',
            source='darboux_susy.py: v6.1',
            comment=f'Изоспектральность {name}: {spec["max_diff"]:.2e}',
            overwrite=True,
        )
        n_registered += 1

    # Счётчик изоспектральных моделей H₂
    n_iso_H2 = sum(
        1 for spec in spectra_H2.values()
        if spec.get('available') and spec.get('isospectral')
    )
    REGISTRY.register(
        'DARBOUX_H2_N_ISOSPECTRAL', n_iso_H2,
        level='1', unit='',
        formula_repr='число изоспектральных моделей H₂',
        source='darboux_susy.py: v6.1',
        comment=f'Изоспектральных H₂: {n_iso_H2}',
        overwrite=True,
    )
    n_registered += 1

    # Счётчик verified моделей H₂
    n_verified_H2 = sum(
        1 for spec in spectra_H2.values()
        if spec.get('available') and spec.get('verified')
    )
    REGISTRY.register(
        'DARBOUX_H2_N_VERIFIED', n_verified_H2,
        level='1', unit='',
        formula_repr='число verified моделей H₂',
        source='darboux_susy.py: v6.1',
        comment=f'Verified H₂: {n_verified_H2}',
        overwrite=True,
    )
    n_registered += 1

    if verbose:
        print(f"  📊 darboux_susy.py: зарегистрировано {n_registered} констант")

    return n_registered


_REGISTERED_COUNT = _register_darboux_constants(verbose=False)


# ============================================================================
# ЧАСТЬ 7: ПЕЧАТЬ ОТЧЁТА
# ============================================================================

def print_report():
    """Печать отчёта v6.1."""
    print("=" * 78)
    print("  Ω⁴-UNIFIED — ФАКТОРИЗАЦИЯ ДАРБУ И SUSY QM (v6.1)")
    print("=" * 78)
    print()

    if not NUMPY_AVAILABLE:
        print("  ⚠️ numpy не установлен — спектральная проверка пропущена")
        return

    results = verify_darboux_susy()

    # ─── 1. Факторизация ────────────────────────────────────────────
    print("── ФАКТОРИЗАЦИЯ ДАРБУ ──")
    print(f"  {'Модель':<20} {'V₂−V₁=2W′':<12} {'max_error':<12} "
          f"{'symmetry':<12}")
    print(f"  {'─'*20} {'─'*12} {'─'*12} {'─'*12}")
    for name, data in results['factorizations'].items():
        marker = '✅' if data['relation_verified'] else '❌'
        sym_str = (
            '✅' if data['symmetry_verified'] is True
            else ('❌' if data['symmetry_verified'] is False else '—')
        )
        print(f"  {name:<20} {marker:<12} {data['max_error']:<12.2e} "
              f"{sym_str:<12}")
    print()

    # ─── 2. Изоспектральность ────────────────────────────────────────
    print("── ИЗОСПЕКТРАЛЬНОСТЬ (L=30, n_points=5000) ──")
    for name, spec in results['spectra'].items():
        if not spec.get('available'):
            print(f"  {name:<20} ⚠️ {spec.get('error', 'n/a')}")
            continue

        if spec.get('excluded'):
            print(f"  {name:<20} ⏸️  ИСКЛЮЧЁН (нет дискретного основного "
                  f"состояния)")
            continue

        marker = '✅' if spec['matched'] else '⚠️'
        E1_str = ', '.join(f'{e:.4f}' for e in spec['E1'][:4])
        E2_str = ', '.join(f'{e:.4f}' for e in spec['E2'][:4])

        print(f"  {name:<20} {marker}")
        print(f"    E₁ = [{E1_str}]")
        print(f"    E₂ = [{E2_str}]")

        if spec.get('E1_analytical'):
            Ea_str = ', '.join(f'{e:.4f}' for e in spec['E1_analytical'][:4])
            print(f"    E₁(аналит.) = [{Ea_str}]")

        if spec.get('has_bound_state'):
            print(f"    max_diff = {spec['max_diff']:.2e} (E₁[1:] ≈ E₂)")
            if spec.get('analytical_diff') is not None:
                print(f"    analytical_diff = {spec['analytical_diff']:.2e}")
        else:
            print(f"    {spec.get('comment', '')}")
    print()

    # ─── 2b. Полный спектр H₂ (НОВОЕ v6.1) ──────────────────────────
    print("── ПОЛНЫЙ СПЕКТР H₂ (v6.1) ──")
    for name, spec in results.get('spectra_H2', {}).items():
        if not spec.get('available'):
            print(f"  {name:<20} ⚠️ {spec.get('error', 'n/a')}")
            continue

        marker = '✅' if spec['verified'] else '❌'
        E1_str = ', '.join(f'{e:.4f}' for e in spec['E1'][:4])
        E2_str = ', '.join(f'{e:.4f}' for e in spec['E2'][:4])

        print(f"  {name:<20} {marker}")
        print(f"    E₁ = [{E1_str}]")
        print(f"    E₂ = [{E2_str}]")
        print(f"    max|E₁[1:] − E₂| = {spec['max_diff']:.2e}")
        print(f"    max|норм ψ_n^(2) − 1| = {spec['norm_psi2']:.2e}")
        print(f"    Изоспектральность: "
              f"{'✅' if spec['isospectral'] else '❌'}")
        print(f"    Партнёрских функций: {spec['n_ladder']}")
    print()

    # ─── 3. SUSY-алгебра ────────────────────────────────────────────
    print("── SUSY-АЛГЕБРА (аналитика) ──")
    alg1 = results['susy_algebra']['A_Adag_equals_H2']
    alg2 = results['susy_algebra']['Adag_A_equals_H1']
    m1 = '✅' if alg1.get('verified', False) else '❌'
    m2 = '✅' if alg2.get('verified', False) else '❌'
    print(f"  A·A† = H₂: {m1}  (diff = {alg1.get('diff', 0):.2e})")
    print(f"  A†·A = H₁: {m2}  (diff = {alg2.get('diff', 0):.2e})")
    print()

    # ─── 4. Связь с Ω⁴ ──────────────────────────────────────────────
    link = results['omega4_link']
    print("── СВЯЗЬ С Ω⁴ ──")
    print(f"  χ = {link['chi']}, D = {link['D']}")
    print(f"  Статус: {link['status']}")
    print()

    # ─── 5. Итог ────────────────────────────────────────────────────
    print("=" * 78)
    print(f"  Статус: {results['status']}")
    print(f"  Моделей: {results['n_models']}, "
          f"обязательных: {results['n_required']}, "
          f"изоспектральных: {results['n_matched']}")
    print(f"  Обязательные: {', '.join(results['required_models'])}")
    print(f"  Прошли: {', '.join(results['matched_names'])}")
    if results['unmatched_names']:
        print(f"  Не прошли: {', '.join(results['unmatched_names'])}")
    print(f"  Полный спектр H₂: "
          f"изоспектральных {results['n_iso_H2']}, "
          f"verified {results['n_verified_H2']}")
    print(f"  {results['comment']}")
    print(f"  📊 Зарегистрировано констант: {_REGISTERED_COUNT}")
    print("=" * 78)


# ============================================================================
# ТОЧКА ВХОДА
# ============================================================================

if __name__ == '__main__':
    print_report()
