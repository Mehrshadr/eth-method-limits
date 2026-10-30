import Mathlib

/-!
# Triangular fooling-set lemma

Core counting step of the ceiling argument: if a set of coordinates of `I − M0` can be
ordered so that no coordinate "reads" an earlier one (entries below the diagonal vanish)
and none reads itself through `M0` (diagonal entries stay nonzero), then the number of
measurements `rank (I − M0)` is at least the size of that set.

`e : n → ι` picks the coordinates, the linear order on `n` is the topological order.
-/

open Matrix

theorem rank_ge_of_triangular_fooling_set
    {ι n K : Type*} [Fintype ι] [Fintype n] [DecidableEq n] [LinearOrder n] [Field K]
    (A : Matrix ι ι K) (e : n → ι)
    (hlow : ∀ i j, j < i → A (e i) (e j) = 0)
    (hdiag : ∀ i, A (e i) (e i) ≠ 0) :
    Fintype.card n ≤ A.rank := by
  set B : Matrix n n K := A.submatrix e e with hB
  have htri : B.IsUpperTriangular := fun i j hij => by simpa [hB] using hlow i j hij
  have hdet : B.det ≠ 0 := by
    rw [det_of_isUpperTriangular htri]
    exact Finset.prod_ne_zero_iff.mpr fun i _ => by simpa [hB] using hdiag i
  have hunit : IsUnit B := (isUnit_iff_isUnit_det B).mpr (Ne.isUnit hdet)
  calc Fintype.card n = B.rank := (rank_of_isUnit B hunit).symm
    _ ≤ A.rank := rank_submatrix_le A e e

/-- The form used in the argument: `A = I − M0`, where `M0` is zero on every pair the
pattern forbids. Coordinates of `S` read only later coordinates of `S` and never themselves. -/
theorem measurements_ge_fooling_set
    {ι n K : Type*} [Fintype ι] [DecidableEq ι] [Fintype n] [DecidableEq n] [LinearOrder n]
    [Field K] (M0 : Matrix ι ι K) (e : n → ι) (he : Function.Injective e)
    (hforward : ∀ i j, j < i → M0 (e i) (e j) = 0)
    (hnoloop : ∀ i, M0 (e i) (e i) = 0) :
    Fintype.card n ≤ (1 - M0).rank := by
  apply rank_ge_of_triangular_fooling_set (1 - M0) e
  · intro i j hij
    have hne : e i ≠ e j := fun h => (ne_of_lt hij).symm (he h)
    simp [one_apply_ne hne, hforward i j hij]
  · intro i
    simp [hnoloop i]

#print axioms rank_ge_of_triangular_fooling_set
#print axioms measurements_ge_fooling_set
