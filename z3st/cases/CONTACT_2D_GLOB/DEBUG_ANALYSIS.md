# Analyse : Pourquoi After_contact a plusieurs points par rayon vs. Cas 9

## Résumé
**Le cas 9 (single cylinder)** capture ~1-2 valeurs uniques de contrainte par rayon. **After_contact** en capturait 33 (avant fix) puis 5 (après optimisation y_tol). La raison est la différence dans la **structure du maillage**.

---

## Comparaison des Maillages

| Paramètre | Cas 9 | After_contact |
|-----------|-------|---------------|
| **Dimensions axiales** | Lz = 0.50 m | Ly = 0.30 m |
| **Divisions axiales (ny)** | 81 | **500** |
| **Densité (éléments/mètre)** | 162 /m | **1667 /m** |
| **Largeur (radiale)** | Ro - Ri = 0.01 m | Ro - Ri = 0.004 m |
| **Divisions radiales (nx)** | 41 | 40 |

---

## Extraction Mid-Plane : Impact de la Tolérance

### Cas 9 (single cylinder)
- Tolérance: **z_tol = 0.01 m** (constant)
- Hauteur totale: 0.50 m
- Couverture: 0.01 / 0.50 = **2% de la hauteur**
- **Éléments capturés par rayon**: ~81 × (0.01/0.50) = **1.6 → arrondi à 4 éléments par rayon unique**
- Points extraits: 160 total → 40 radii uniques = **4 points/rayon**
- ✅ **Résultat**: PASS (erreur L2 = 1.67e-04)

### After_contact (avant optimisation)
- Tolérance: **y_tol = 0.01 m** (constant)
- Hauteur totale: 0.30 m
- Couverture: 0.01 / 0.30 = **3.3% de la hauteur**
- **Éléments capturés par rayon**: ~500 × (0.01/0.30) = **16.7 → 33 points/rayon**
- Points extraits: 2574 total → 83 radii uniques = **31 points/rayon**
- ❌ **Problème**: 31 valeurs de contrainte pour 1 seul rayon!

### After_contact (après optimisation)
- Tolérance: **y_tol = 0.0016 m** (scaled: 0.01 × 81/500)
- Couverture: 0.0016 / 0.30 = **0.53% de la hauteur**
- **Éléments capturés par rayon**: ~500 × (0.0016/0.30) = **2.67 → ~5 points/rayon**
- Points extraits: 390 total → 82 radii uniques
- Points restants après agrégation: **82 radii uniques** ✅

---

## Solution Implémentée

### Changement 1: Tolérance Adaptée
```python
# Avant
y_tol = 0.01  # Capturait trop d'éléments

# Après
y_tol = 0.0016  # Scaled: 0.01 * (81/500) = tolérance équivalente au cas 9
```

### Changement 2: Agrégation par Rayon Unique
```python
def aggregate_by_radius(r_vals, field_vals_3x3):
    """Moyenne les composantes stress/strain par rayon unique."""
    unique_r, inverse_indices = np.unique(r_vals, return_inverse=True)
    n_unique = len(unique_r)
    n_components = field_vals_3x3.shape[1]
    aggregated_field = np.zeros((n_unique, n_components))
    
    for i, r_val in enumerate(unique_r):
        mask = inverse_indices == i
        aggregated_field[i, :] = np.mean(field_vals_3x3[mask, :], axis=0)
    
    return unique_r, aggregated_field
```

**Résultat**: 
- Avant agrégation: 390 points avec doublons
- Après agrégation: 82 radii uniques (outer), 100 (inner) = **1 valeur par rayon** ✅

---

## Pourquoi le Cas 9 Fonctionne Malgré les Doublons

Le cas 9 a aussi 4 points/rayon mais **passe les tests**. Raison:
1. Les 4 points par rayon sont des **doublons exacts** (même rayon x = constant)
2. La solution analytique Lamé donne la même valeur pour chaque rayon répété: $\sigma_{rr}(r) = A - B/r^2$
3. Les tableaux restent **cohérents** (4 points, 4 valeurs analytiques identiques)
4. L'erreur L2 moyenne = très faible

**After_contact (avant fix)** avait 31 doublons par rayon, ce qui créait une **variabilité indésirable** dans les profils de contrainte tracés.

---

## Fichiers Modifiés

1. **[After_contact/non-regression.py](../../After_contact/non-regression.py)**
   - Ligne ~62: Réduit `y_tol` de 0.01 à 0.0016
   - Lignes ~92-115: Ajoute la fonction `aggregate_by_radius()`
   - Lignes ~120-145: Applique agrégation au stress et strain

2. **[9_thick_cylindrical_shell_GPS_2D/non-regression.py](../../9_thick_cylindrical_shell_GPS_2D/non-regression.py)**
   - Lignes ~62-77: Ajoute debug pour visualiser les doublons
   - Résultat: 40 radii uniques, 4 points/rayon (cohérent)

---

## Conclusion

| Aspect | Cas 9 | After_contact Before | After_contact After |
|--------|-------|---------------------|---------------------|
| **Points mid-plane** | 160 | 2574 | 390 |
| **Radii uniques** | 40 | 83 | 82 |
| **Points/rayon** | 4 | 31 | 5 (→1 après agrégation) |
| **Cohérence** | ✅ Bonne | ❌ Mauvaise | ✅ Bonne |
| **Non-regression** | ✅ PASS | ❌ FAIL | ⚠️ FAIL (mais structure OK) |

**Observation clé**: La structure du maillage, pas les doublons per se, cause la différence. Avec 500 divisions axiales (vs 81), une tolérance fixe capture beaucoup plus d'éléments. La solution: **adapter la tolérance au rapport de densité de maillage + agréger par rayon unique**.
