#!/usr/bin/env python3
"""
Theoretical Fragment Ion Calculator for O-GlcNAc Glycopeptides

This module calculates theoretical m/z values for peptide fragment ions
including b, y, c, z ions with modifications, neutral losses, and
glycan-specific ions for EThcD spectrum annotation.

Author: Claude Code Assistant
Date: 2026-01-21
"""

import numpy as np
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional, Set
import json
import re

# =============================================================================
# CONSTANTS: Monoisotopic Masses
# =============================================================================

# Amino acid residue masses (monoisotopic)
AA_MASSES = {
    'A': 71.03711,   # Alanine
    'R': 156.10111,  # Arginine
    'N': 114.04293,  # Asparagine
    'D': 115.02694,  # Aspartic acid
    'C': 103.00919,  # Cysteine
    'E': 129.04259,  # Glutamic acid
    'Q': 128.05858,  # Glutamine
    'G': 57.02146,   # Glycine
    'H': 137.05891,  # Histidine
    'I': 113.08406,  # Isoleucine
    'L': 113.08406,  # Leucine
    'K': 128.09496,  # Lysine
    'M': 131.04049,  # Methionine
    'F': 147.06841,  # Phenylalanine
    'P': 97.05276,   # Proline
    'S': 87.03203,   # Serine
    'T': 101.04768,  # Threonine
    'W': 186.07931,  # Tryptophan
    'Y': 163.06333,  # Tyrosine
    'V': 99.06841,   # Valine
}

# Modification masses
MOD_MASSES = {
    'TMT6plex': 229.1629,           # TMT6plex tag
    'HexNAc': 203.0794,             # N-acetylhexosamine (O-GlcNAc)
    'HexNAc_TMT': 528.2859,         # O-GlcNAc with TMT reporter
    'Carbamidomethyl': 57.02146,    # Cysteine alkylation
    'Oxidation': 15.9949,           # Methionine oxidation
}

# Common masses
PROTON = 1.007276
H2O = 18.010565
NH3 = 17.026549
CO = 27.994915
ELECTRON = 0.000549

# Neutral loss masses
NEUTRAL_LOSSES = {
    'H2O': 18.010565,
    'NH3': 17.026549,
    'HexNAc_TMT': 528.2859,  # Loss of glycan with TMT
    'HexNAc': 300.1308,       # Loss of glycan (oxonium form)
}

# Oxonium ions (glycan diagnostic ions) - B-type ions
OXONIUM_IONS = {
    'HexNAc_TMT': 529.2937,      # HexNAc with TMT oxonium
    'HexNAc+': 300.1308,         # HexNAc + H2O + H
    'HexNAc': 204.0864,          # HexNAc oxonium
    'HexNAc-H2O': 186.0760,      # HexNAc - H2O
    'HexNAc-2H2O': 168.0652,     # HexNAc - 2H2O
    'HexNAc_frag1': 144.0652,    # HexNAc fragment
    'HexNAc_frag2': 138.0546,    # HexNAc fragment
}

# Glycan Y ion definitions (peptide + partial glycan)
# For O-GlcNAc with TMT: Y0 = peptide only, Y1 = peptide + HexNAc+TMT (intact)
# Mass values represent the LOSS from intact glycopeptide to get each Y ion
GLYCAN_Y_LOSSES = {
    'HexNAc_TMT': {
        # O-GlcNAc + TMT (528.2859 Da total)
        'Y0': 528.2859,           # Complete loss: peptide only
        'Y0-H2O': 528.2859 - 18.0106,  # Y0 with water loss from peptide
        'Y*': 203.0794,           # Loss of HexNAc only, TMT stays on peptide (rare)
    },
    'HexNAc': {
        # O-GlcNAc without TMT (203.0794 Da)
        'Y0': 203.0794,           # Complete loss: peptide only
        'Y0-H2O': 203.0794 - 18.0106,  # Y0 with water loss
    },
}

# Ion type mass adjustments (added to residue sum)
# For a fragment containing residues with total mass M:
ION_ADJUSTMENTS = {
    'b': PROTON,                           # +H+
    'y': H2O + PROTON,                     # +H2O +H+
    'c': NH3 + PROTON,                     # +NH3 +H+ (or +NH2 to peptide)
    'z': -NH3 + H2O + 2 * PROTON,          # z+H (even-electron, H rearrangement)
}

# =============================================================================
# DATA CLASSES
# =============================================================================

@dataclass
class Modification:
    """Represents a modification on a peptide."""
    position: int       # 0 = N-term, -1 = C-term, 1-based for residues
    residue: str        # 'N-term', 'C-term', or amino acid letter
    mass: float         # Modification mass
    name: str = ""      # Optional name


@dataclass
class TheoreticalIon:
    """Represents a theoretical fragment ion."""
    ion_type: str       # 'b', 'y', 'c', 'z', 'Y', 'oxonium', 'precursor'
    ion_number: int     # Ion number (e.g., b3 = 3)
    charge: int         # Charge state
    mz: float           # Theoretical m/z
    sequence: str       # Fragment sequence
    neutral_loss: str = ""  # Neutral loss type if any
    annotation: str = ""    # Full annotation string


@dataclass
class MatchedIon(TheoreticalIon):
    """Represents a matched ion with experimental data."""
    exp_mz: float = 0.0
    exp_intensity: float = 0.0
    mass_error_ppm: float = 0.0


# =============================================================================
# FRAGMENT CALCULATOR CLASS
# =============================================================================

class FragmentCalculator:
    """
    Calculates theoretical fragment ions for glycopeptides.

    Supports:
    - b, y ions (HCD)
    - c, z ions (ETD)
    - Y ions (precursor with glycan)
    - Charge-reduced precursor species
    - Oxonium ions
    - Neutral losses (H2O, NH3, glycan)
    """

    def __init__(self,
                 peptide: str,
                 modifications: List[Dict],
                 precursor_charge: int,
                 max_fragment_charge: int = 2):
        """
        Initialize the fragment calculator.

        Args:
            peptide: Peptide sequence (single letter amino acids)
            modifications: List of dicts with 'position', 'residue', 'mass'
            precursor_charge: Precursor ion charge state
            max_fragment_charge: Maximum fragment ion charge (default 2)
        """
        self.peptide = peptide.upper()
        self.length = len(peptide)
        self.precursor_charge = precursor_charge
        self.max_fragment_charge = min(max_fragment_charge, precursor_charge - 1)

        # Parse modifications
        self.modifications = []
        self.mod_by_position = {}  # position -> Modification

        for mod in modifications:
            m = Modification(
                position=mod['position'],
                residue=mod.get('residue', ''),
                mass=mod['mass'],
                name=mod.get('name', '')
            )
            self.modifications.append(m)
            self.mod_by_position[m.position] = m

        # Calculate residue masses including modifications
        self._calculate_residue_masses()

        # Calculate precursor mass
        self.precursor_mass = self._calculate_precursor_mass()
        self.precursor_mz = (self.precursor_mass + precursor_charge * PROTON) / precursor_charge

    def _calculate_residue_masses(self):
        """Calculate mass of each residue including any modifications."""
        self.residue_masses = []

        for i, aa in enumerate(self.peptide):
            pos = i + 1  # 1-based position
            mass = AA_MASSES.get(aa, 0)

            # Add modification mass if present
            if pos in self.mod_by_position:
                mass += self.mod_by_position[pos].mass

            self.residue_masses.append(mass)

        # Store N-term and C-term modifications separately
        self.nterm_mod_mass = self.mod_by_position.get(0, Modification(0, '', 0)).mass
        self.cterm_mod_mass = self.mod_by_position.get(-1, Modification(-1, '', 0)).mass

    def _calculate_precursor_mass(self) -> float:
        """Calculate the neutral monoisotopic mass of the precursor."""
        # Sum of residues + H2O (terminal groups) + modifications
        mass = sum(self.residue_masses) + H2O
        mass += self.nterm_mod_mass
        mass += self.cterm_mod_mass
        return mass

    def _get_glycan_position(self) -> Optional[int]:
        """Find the position of the glycan modification (528.2859 mass)."""
        for pos, mod in self.mod_by_position.items():
            if abs(mod.mass - MOD_MASSES['HexNAc_TMT']) < 0.01:
                return pos
        return None

    def calculate_b_ions(self, charges: List[int] = None) -> List[TheoreticalIon]:
        """
        Calculate b ions (N-terminal, HCD).
        b_n = sum of first n residues + N-term mod + H+
        """
        if charges is None:
            charges = list(range(1, self.max_fragment_charge + 1))

        ions = []
        cumulative_mass = self.nterm_mod_mass

        for i in range(self.length - 1):  # b1 to b(n-1)
            cumulative_mass += self.residue_masses[i]
            ion_mass = cumulative_mass + ION_ADJUSTMENTS['b'] - PROTON  # Neutral mass

            for charge in charges:
                mz = (ion_mass + charge * PROTON) / charge
                ion = TheoreticalIon(
                    ion_type='b',
                    ion_number=i + 1,
                    charge=charge,
                    mz=mz,
                    sequence=self.peptide[:i+1],
                    annotation=f"b{i+1}{'⁺' * charge}"
                )
                ions.append(ion)

        return ions

    def calculate_y_ions(self, charges: List[int] = None) -> List[TheoreticalIon]:
        """
        Calculate y ions (C-terminal, HCD).
        y_n = sum of last n residues + C-term mod + H2O + H+
        """
        if charges is None:
            charges = list(range(1, self.max_fragment_charge + 1))

        ions = []
        cumulative_mass = self.cterm_mod_mass

        for i in range(self.length - 1):  # y1 to y(n-1)
            cumulative_mass += self.residue_masses[self.length - 1 - i]
            ion_mass = cumulative_mass + ION_ADJUSTMENTS['y'] - PROTON  # Neutral mass

            for charge in charges:
                mz = (ion_mass + charge * PROTON) / charge
                ion = TheoreticalIon(
                    ion_type='y',
                    ion_number=i + 1,
                    charge=charge,
                    mz=mz,
                    sequence=self.peptide[-(i+1):],
                    annotation=f"y{i+1}{'⁺' * charge}"
                )
                ions.append(ion)

        return ions

    def calculate_c_ions(self, charges: List[int] = None) -> List[TheoreticalIon]:
        """
        Calculate c ions (N-terminal, ETD).
        c_n = b_n + NH3 = sum of first n residues + N-term mod + NH3 + H+
        """
        if charges is None:
            charges = list(range(1, self.max_fragment_charge + 1))

        ions = []
        cumulative_mass = self.nterm_mod_mass

        for i in range(self.length - 1):  # c1 to c(n-1)
            cumulative_mass += self.residue_masses[i]
            ion_mass = cumulative_mass + ION_ADJUSTMENTS['c'] - PROTON  # Neutral mass

            for charge in charges:
                mz = (ion_mass + charge * PROTON) / charge
                ion = TheoreticalIon(
                    ion_type='c',
                    ion_number=i + 1,
                    charge=charge,
                    mz=mz,
                    sequence=self.peptide[:i+1],
                    annotation=f"c{i+1}{'⁺' * charge}"
                )
                ions.append(ion)

        return ions

    def calculate_z_ions(self, charges: List[int] = None) -> List[TheoreticalIon]:
        """
        Calculate z ions (C-terminal, ETD).
        z_n = y_n - NH3 (z-radical ion)
        """
        if charges is None:
            charges = list(range(1, self.max_fragment_charge + 1))

        ions = []
        cumulative_mass = self.cterm_mod_mass

        for i in range(self.length - 1):  # z1 to z(n-1)
            cumulative_mass += self.residue_masses[self.length - 1 - i]
            ion_mass = cumulative_mass + ION_ADJUSTMENTS['z'] - PROTON  # Neutral mass

            for charge in charges:
                mz = (ion_mass + charge * PROTON) / charge
                ion = TheoreticalIon(
                    ion_type='z',
                    ion_number=i + 1,
                    charge=charge,
                    mz=mz,
                    sequence=self.peptide[-(i+1):],
                    annotation=f"z{i+1}{'⁺' * charge}"
                )
                ions.append(ion)

        return ions

    def calculate_Y_ions(self, charges: List[int] = None) -> List[TheoreticalIon]:
        """
        Calculate glycan Y ions (peptide + partial/no glycan).

        Y ions in glycopeptide MS represent the peptide backbone with
        varying amounts of glycan attached:
        - Y0: Peptide only (complete glycan loss)
        - Y1: Peptide + 1 monosaccharide (for complex glycans)
        - Y(intact): Full glycopeptide

        For O-GlcNAc (single monosaccharide):
        - Y0 = peptide backbone (glycan completely lost)
        - Y1 = intact glycopeptide
        """
        if charges is None:
            charges = list(range(1, self.precursor_charge + 1))

        ions = []

        # Determine glycan type and get Y ion losses
        glycan_type = self._get_glycan_type()
        glycan_position = self._get_glycan_position()

        if glycan_type and glycan_type in GLYCAN_Y_LOSSES:
            y_losses = GLYCAN_Y_LOSSES[glycan_type]

            # Calculate each Y ion type
            for y_name, loss_mass in y_losses.items():
                y_mass = self.precursor_mass - loss_mass

                for charge in charges:
                    mz = (y_mass + charge * PROTON) / charge

                    # Determine Y ion number (0 = no glycan, 1 = full glycan for O-GlcNAc)
                    if 'Y0' in y_name:
                        ion_number = 0
                        annotation = f"{y_name} {charge}+"
                    else:
                        ion_number = 1
                        annotation = f"Y* {charge}+"

                    ion = TheoreticalIon(
                        ion_type='Y',
                        ion_number=ion_number,
                        charge=charge,
                        mz=mz,
                        sequence=self.peptide,
                        annotation=annotation
                    )
                    ions.append(ion)

        # Always add intact glycopeptide (Y1 or Y-full) at different charge states
        for charge in charges:
            mz = (self.precursor_mass + charge * PROTON) / charge
            ion = TheoreticalIon(
                ion_type='Y',
                ion_number=1 if glycan_type else 0,
                charge=charge,
                mz=mz,
                sequence=self.peptide,
                annotation=f"Y1 {charge}+" if glycan_type else f"[M+{charge}H]{charge}+"
            )
            ions.append(ion)

        return ions

    def _get_glycan_type(self) -> Optional[str]:
        """Determine the type of glycan modification."""
        for pos, mod in self.mod_by_position.items():
            if pos == 0 or pos == -1:
                continue  # skip N/C-term mods
            # HexNAc + TMT = 528.2859 Da
            if abs(mod.mass - MOD_MASSES['HexNAc_TMT']) < 0.1:
                return 'HexNAc_TMT'
            # HexNAc only (metabolic labeling) = 203.0794 Da
            elif abs(mod.mass - MOD_MASSES['HexNAc']) < 0.1:
                return 'HexNAc'
            # OPair-reported O-GlcNAc = 299.1230 Da (HexNAc_TMT - TMT)
            elif abs(mod.mass - 299.1230) < 0.1:
                return 'HexNAc'
        return None

    def calculate_charge_reduced_precursor(self) -> List[TheoreticalIon]:
        """
        Calculate charge-reduced precursor species from ETD.

        These are the intact glycopeptide at reduced charge states (ETnoD).
        Labeled as Y1 ions in glycoproteomics nomenclature since the glycan
        is still attached: Y1 = peptide + glycan (intact glycopeptide).
        """
        ions = []
        glycan_type = self._get_glycan_type()

        for reduced_charge in range(1, self.precursor_charge):
            mz = (self.precursor_mass + reduced_charge * PROTON) / reduced_charge
            # Use Y1 label for glycopeptides (glycan intact), generic label otherwise
            if glycan_type:
                annotation = f"Y1 {reduced_charge}+"
                ion_type = 'Y'
                ion_number = 1
            else:
                annotation = f"[M+{reduced_charge}H]{reduced_charge}+"
                ion_type = 'precursor'
                ion_number = 0
            ion = TheoreticalIon(
                ion_type=ion_type,
                ion_number=ion_number,
                charge=reduced_charge,
                mz=mz,
                sequence=self.peptide,
                annotation=annotation
            )
            ions.append(ion)

        return ions

    def calculate_precursor_isotopes(self, n_isotopes: int = 4) -> List[TheoreticalIon]:
        """Calculate precursor isotope peaks."""
        ions = []
        isotope_spacing = 1.003355 / self.precursor_charge

        for i in range(n_isotopes):
            mz = self.precursor_mz + i * isotope_spacing
            ion = TheoreticalIon(
                ion_type='precursor',
                ion_number=i,
                charge=self.precursor_charge,
                mz=mz,
                sequence=self.peptide,
                annotation=f"[M+{self.precursor_charge}H]{self.precursor_charge}+ iso{i}"
            )
            ions.append(ion)

        return ions

    def calculate_oxonium_ions(self) -> List[TheoreticalIon]:
        """Calculate glycan oxonium (diagnostic) ions."""
        ions = []

        for name, mz in OXONIUM_IONS.items():
            ion = TheoreticalIon(
                ion_type='oxonium',
                ion_number=0,
                charge=1,
                mz=mz,
                sequence='',
                annotation=name
            )
            ions.append(ion)

        return ions

    def calculate_neutral_loss_ions(self,
                                    base_ions: List[TheoreticalIon],
                                    loss_types: List[str] = None) -> List[TheoreticalIon]:
        """
        Calculate neutral loss ions from base fragment ions.

        Args:
            base_ions: List of base fragment ions
            loss_types: List of neutral loss types ('H2O', 'NH3', 'HexNAc_TMT', 'HexNAc')
        """
        if loss_types is None:
            loss_types = ['H2O', 'NH3']

        glycan_pos = self._get_glycan_position()
        ions = []

        for base_ion in base_ions:
            for loss_type in loss_types:
                loss_mass = NEUTRAL_LOSSES.get(loss_type, 0)
                if loss_mass == 0:
                    continue

                # For glycan loss, only apply if the fragment contains the glycan
                if loss_type in ['HexNAc_TMT', 'HexNAc'] and glycan_pos:
                    # Check if this fragment contains the glycan position
                    if base_ion.ion_type in ['b', 'c']:
                        if base_ion.ion_number < glycan_pos:
                            continue  # Fragment doesn't include glycan
                    elif base_ion.ion_type in ['y', 'z']:
                        if base_ion.ion_number < (self.length - glycan_pos + 1):
                            continue  # Fragment doesn't include glycan

                # Calculate neutral loss m/z
                new_mz = base_ion.mz - loss_mass / base_ion.charge

                if new_mz > 0:
                    ion = TheoreticalIon(
                        ion_type=base_ion.ion_type,
                        ion_number=base_ion.ion_number,
                        charge=base_ion.charge,
                        mz=new_mz,
                        sequence=base_ion.sequence,
                        neutral_loss=loss_type,
                        annotation=f"{base_ion.annotation}-{loss_type}"
                    )
                    ions.append(ion)

        return ions

    def calculate_all_ions(self,
                          include_neutral_losses: bool = True,
                          neutral_loss_types: List[str] = None) -> Dict[str, List[TheoreticalIon]]:
        """
        Calculate all theoretical ions for the peptide.

        Returns dict organized by ion type.
        """
        if neutral_loss_types is None:
            neutral_loss_types = ['H2O', 'NH3', 'HexNAc_TMT', 'HexNAc']

        # Charge-reduced precursors are Y1 ions for glycopeptides
        y_ions = self.calculate_Y_ions()
        cr_ions = self.calculate_charge_reduced_precursor()

        result = {
            'b': self.calculate_b_ions(),
            'y': self.calculate_y_ions(),
            'c': self.calculate_c_ions(),
            'z': self.calculate_z_ions(),
            'Y': y_ions + cr_ions,
            'precursor': self.calculate_precursor_isotopes(),
            'oxonium': self.calculate_oxonium_ions(),
        }

        if include_neutral_losses:
            # Add neutral losses for backbone ions
            for ion_type in ['b', 'y', 'c', 'z']:
                nl_ions = self.calculate_neutral_loss_ions(result[ion_type], neutral_loss_types)
                result[f'{ion_type}_NL'] = nl_ions

        return result

    def get_all_ions_flat(self, **kwargs) -> List[TheoreticalIon]:
        """Get all ions as a flat list."""
        all_ions = self.calculate_all_ions(**kwargs)
        flat_list = []
        for ion_list in all_ions.values():
            flat_list.extend(ion_list)
        return flat_list


# =============================================================================
# DEISOTOPING
# =============================================================================

ISOTOPE_SPACING = 1.003355  # C13-C12 mass difference


def deisotope(exp_mz: np.ndarray,
              exp_intensity: np.ndarray,
              max_charge: int = 4,
              tolerance_ppm: float = 15.0,
              min_isotope_ratio: float = 0.01,
              max_isotope_ratio: float = 2.0,
              ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Remove isotope peaks from a centroided spectrum, keeping monoisotopic peaks.

    For each peak, checks whether a peak exists ~1.003/z Da higher (M+1) that
    could be an isotope. If found, marks the higher peak as an isotope and
    removes it.

    Args:
        exp_mz: Experimental m/z array (sorted ascending)
        exp_intensity: Experimental intensity array
        max_charge: Maximum charge state to consider for isotope spacing
        tolerance_ppm: Mass tolerance for isotope peak matching
        min_isotope_ratio: Minimum M+1/M+0 intensity ratio to be an isotope
        max_isotope_ratio: Maximum M+1/M+0 intensity ratio (for large
            peptides where M+1 can exceed M+0)

    Returns:
        (deisotoped_mz, deisotoped_intensity) with isotope peaks removed
    """
    if len(exp_mz) == 0:
        return exp_mz, exp_intensity

    is_isotope = np.zeros(len(exp_mz), dtype=bool)

    for i in range(len(exp_mz)):
        if is_isotope[i]:
            continue
        for z in range(1, max_charge + 1):
            spacing = ISOTOPE_SPACING / z
            target_mz = exp_mz[i] + spacing
            tol = target_mz * tolerance_ppm / 1e6
            # Search forward for M+1 peak
            for j in range(i + 1, len(exp_mz)):
                if exp_mz[j] > target_mz + tol:
                    break
                if abs(exp_mz[j] - target_mz) <= tol and not is_isotope[j]:
                    ratio = exp_intensity[j] / exp_intensity[i] if exp_intensity[i] > 0 else 0
                    if min_isotope_ratio <= ratio <= max_isotope_ratio:
                        is_isotope[j] = True
                        # Also check M+2
                        target_m2 = exp_mz[i] + 2 * spacing
                        tol_m2 = target_m2 * tolerance_ppm / 1e6
                        for k in range(j + 1, len(exp_mz)):
                            if exp_mz[k] > target_m2 + tol_m2:
                                break
                            if abs(exp_mz[k] - target_m2) <= tol_m2 and not is_isotope[k]:
                                ratio2 = exp_intensity[k] / exp_intensity[i] if exp_intensity[i] > 0 else 0
                                if ratio2 <= max_isotope_ratio:
                                    is_isotope[k] = True
                                    break
                    break

    keep = ~is_isotope
    return exp_mz[keep], exp_intensity[keep]


# =============================================================================
# CHARGE-REDUCED PRECURSOR EXCLUSION (EThcD)
# =============================================================================

GLYCAN_RESIDUE_MASSES = {
    'HexNAc': 203.07937,
    'Hex': 162.05282,
    'Fuc': 146.05791,
    'NeuAc': 291.09542,
}


def build_charge_reduced_exclusions(
    precursor_neutral_mass: float,
    charge: int,
    glycan_mass: float,
    n_isotopes: int = 4,
) -> List[float]:
    """
    Build list of m/z values for charge-reduced precursor species (ETnoD)
    and their common neutral loss satellites.

    In ETD/EThcD, the precursor can capture an electron without backbone
    fragmentation, producing charge-reduced species at z-1, z-2, etc.
    These can be falsely matched as c/z fragment ions.

    Returns list of m/z values to exclude from c/z matching.
    """
    exclusions = []
    neutral_losses = [
        0.0,
        203.07937,   # -HexNAc
        162.05282,   # -Hex
        291.09542,   # -NeuAc
    ]

    for reduced_z in range(1, charge):
        for nl_mass in neutral_losses:
            if nl_mass > glycan_mass + 0.01:
                continue
            for iso in range(n_isotopes):
                cr_mz = (precursor_neutral_mass - nl_mass
                         + reduced_z * PROTON
                         + iso * ISOTOPE_SPACING) / reduced_z
                exclusions.append(cr_mz)

    return exclusions


def filter_charge_reduced(
    matched_ions: List['MatchedIon'],
    exclusion_mzs: List[float],
    exclusion_ppm: float = 15.0,
) -> List['MatchedIon']:
    """
    Remove matched c/z ions that overlap with charge-reduced precursor species.

    Args:
        matched_ions: List of matched ions from match_peaks
        exclusion_mzs: List of m/z values to exclude (from build_charge_reduced_exclusions)
        exclusion_ppm: Tolerance in ppm for exclusion matching

    Returns:
        Filtered list with charge-reduced artifacts removed
    """
    if not exclusion_mzs:
        return matched_ions

    exclusion_arr = np.array(exclusion_mzs)
    filtered = []
    for m in matched_ions:
        if m.ion_type in ('c', 'z'):
            ppm_diffs = np.abs(m.exp_mz - exclusion_arr) / exclusion_arr * 1e6
            if np.min(ppm_diffs) < exclusion_ppm:
                continue
        filtered.append(m)
    return filtered


# =============================================================================
# S/N CALCULATION
# =============================================================================

def load_noise_cache(noise_mzml_path: str) -> Dict[int, Dict[str, np.ndarray]]:
    """
    Load noise profiles from ThermoRawFileParser -N mzML file.

    Returns dict keyed by scan number with 'noise_mz' and 'noise_int' arrays.
    """
    from pyteomics import mzml
    cache = {}
    reader = mzml.MzML(str(noise_mzml_path))
    for spec in reader:
        scan_id = spec.get('id', '')
        if 'scan=' in scan_id:
            scan_num = int(scan_id.split('scan=')[-1])
        else:
            continue
        if 'sampled noise m/z array' in spec:
            cache[scan_num] = {
                'noise_mz': spec['sampled noise m/z array'],
                'noise_int': spec['sampled noise intensity array'],
            }
    reader.close()
    return cache


def get_sn_array(exp_mz: np.ndarray,
                 exp_intensity: np.ndarray,
                 noise_cache: Optional[Dict] = None,
                 scan_num: Optional[int] = None,
                 ) -> np.ndarray:
    """
    Calculate S/N for each peak. Uses instrument noise if available,
    otherwise falls back to median-based estimation.

    Args:
        exp_mz: Experimental m/z array
        exp_intensity: Experimental intensity array
        noise_cache: Dict from load_noise_cache (optional)
        scan_num: Scan number to look up in noise_cache (optional)

    Returns:
        Array of S/N values, one per peak
    """
    if len(exp_mz) == 0:
        return np.array([])

    # Try instrument noise first
    if noise_cache is not None and scan_num is not None:
        noise_data = noise_cache.get(scan_num)
        if noise_data is not None:
            noise_at_peaks = np.interp(exp_mz, noise_data['noise_mz'],
                                        noise_data['noise_int'])
            noise_at_peaks = np.maximum(noise_at_peaks, 1.0)  # avoid div by zero
            return exp_intensity / noise_at_peaks

    # Fallback: median-based noise estimation
    noise_level = max(np.median(exp_intensity), 1.0)
    return exp_intensity / noise_level


# =============================================================================
# PEAK MATCHING
# =============================================================================

def match_peaks(theoretical_ions: List[TheoreticalIon],
                exp_mz: np.ndarray,
                exp_intensity: np.ndarray,
                tolerance_ppm: float = 20.0,
                match_isotopes: bool = True,
                max_isotope: int = 2) -> List[MatchedIon]:
    """
    Match experimental peaks to theoretical ions, including isotope peaks.

    Args:
        theoretical_ions: List of theoretical ions
        exp_mz: Experimental m/z array
        exp_intensity: Experimental intensity array
        tolerance_ppm: Mass tolerance in ppm
        match_isotopes: Whether to also match M+1, M+2 isotope peaks
        max_isotope: Maximum isotope offset to check (default 2 for M+0, M+1, M+2)

    Returns:
        List of matched ions with experimental data
    """
    matched = []
    used_peaks = set()  # Track which peaks have been matched

    # Isotope mass spacing (C13 - C12)
    ISOTOPE_SPACING = 1.003355

    # Sort theoretical ions by priority:
    # 1. Y ions (glycan) first - these are important diagnostic ions
    # 2. Base ions before neutral loss ions
    # 3. Then by m/z
    def ion_sort_key(ion):
        # Priority: Y ions get highest priority (0), then base ions (1), then NL ions (2)
        if ion.ion_type == 'Y':
            type_priority = 0
        elif not ion.neutral_loss:
            type_priority = 1
        else:
            type_priority = 2
        return (type_priority, ion.mz)

    sorted_ions = sorted(theoretical_ions, key=ion_sort_key)

    for ion in sorted_ions:
        theo_mz = ion.mz
        found_match = False

        # Try monoisotopic first, then isotopes if enabled
        isotope_offsets = [0]
        if match_isotopes:
            isotope_offsets = list(range(max_isotope + 1))  # [0, 1, 2]

        for iso_offset in isotope_offsets:
            if found_match:
                break

            # Calculate m/z with isotope offset (divided by charge)
            target_mz = theo_mz + (iso_offset * ISOTOPE_SPACING / ion.charge)
            tol = target_mz * tolerance_ppm / 1e6

            # Find peaks within tolerance
            matches = np.where(np.abs(exp_mz - target_mz) <= tol)[0]

            if len(matches) > 0:
                # Find the closest match that hasn't been used
                for idx in matches[np.argsort(np.abs(exp_mz[matches] - target_mz))]:
                    if idx not in used_peaks:
                        # Calculate error relative to theoretical monoisotopic
                        error_ppm = (exp_mz[idx] - target_mz) / target_mz * 1e6

                        # Update annotation if isotope match
                        annotation = ion.annotation
                        if iso_offset > 0:
                            annotation = f"{ion.annotation}+{iso_offset}"

                        matched_ion = MatchedIon(
                            ion_type=ion.ion_type,
                            ion_number=ion.ion_number,
                            charge=ion.charge,
                            mz=ion.mz,  # Keep theoretical monoisotopic m/z
                            sequence=ion.sequence,
                            neutral_loss=ion.neutral_loss,
                            annotation=annotation,
                            exp_mz=exp_mz[idx],
                            exp_intensity=exp_intensity[idx],
                            mass_error_ppm=error_ppm
                        )
                        matched.append(matched_ion)
                        used_peaks.add(idx)
                        found_match = True
                        break

    return matched


# =============================================================================
# FALSE MATCH RATE CALCULATION
# =============================================================================

@dataclass
class FalseMatchRate:
    """Results of false match rate estimation."""
    fmr_peaks: float          # False match rate based on peak count
    fmr_intensity: float      # False match rate based on intensity
    matched_peaks: int        # Number of matched peaks in true spectrum
    matched_intensity: float  # Sum of matched intensity in true spectrum
    avg_random_peaks: float   # Average matched peaks in shifted spectra
    avg_random_intensity: float  # Average matched intensity in shifted spectra
    n_shifts: int             # Number of spectrum shifts used


def calculate_false_match_rate(
    theoretical_ions: List[TheoreticalIon],
    exp_mz: np.ndarray,
    exp_intensity: np.ndarray,
    tolerance_ppm: float = 20.0,
    shift_range: float = 25.0,
    shift_step: float = 1.0
) -> FalseMatchRate:
    """
    Calculate false match rate using spectrum shifting method.

    This implements the method from Schulte et al. (Anal. Chem. 2025):
    The spectrum is shifted by π ± shift_range Th in shift_step increments.
    The π offset prevents false matches from isotope patterns.

    Args:
        theoretical_ions: List of theoretical fragment ions
        exp_mz: Experimental m/z array
        exp_intensity: Experimental intensity array
        tolerance_ppm: Mass tolerance for matching (default 20 ppm)
        shift_range: Range of shifts in Th (default ±25 Th)
        shift_step: Step size for shifts (default 1 Th)

    Returns:
        FalseMatchRate object with peak and intensity-based FMR
    """
    # First, match the true (unshifted) spectrum
    true_matched = match_peaks(
        theoretical_ions, exp_mz, exp_intensity,
        tolerance_ppm, match_isotopes=False  # Disable isotope matching for FMR
    )

    true_peak_count = len(true_matched)
    true_intensity_sum = sum(ion.exp_intensity for ion in true_matched)

    # If no matches in true spectrum, return zero FMR
    if true_peak_count == 0:
        return FalseMatchRate(
            fmr_peaks=0.0,
            fmr_intensity=0.0,
            matched_peaks=0,
            matched_intensity=0.0,
            avg_random_peaks=0.0,
            avg_random_intensity=0.0,
            n_shifts=0
        )

    # Generate shift offsets: π ± shift_range in shift_step increments
    # π is added to prevent false matches from isotope patterns (spacing ~1 Da)
    pi_offset = np.pi
    shifts = np.arange(-shift_range, shift_range + shift_step, shift_step) + pi_offset

    # Track matches for each shifted spectrum
    shifted_peak_counts = []
    shifted_intensity_sums = []

    for shift in shifts:
        # Shift the experimental m/z values
        shifted_mz = exp_mz + shift

        # Match against theoretical ions
        shifted_matched = match_peaks(
            theoretical_ions, shifted_mz, exp_intensity,
            tolerance_ppm, match_isotopes=False
        )

        shifted_peak_counts.append(len(shifted_matched))
        shifted_intensity_sums.append(sum(ion.exp_intensity for ion in shifted_matched))

    # Calculate average matches across all shifts
    avg_random_peaks = np.mean(shifted_peak_counts)
    avg_random_intensity = np.mean(shifted_intensity_sums)

    # Calculate false match rates
    fmr_peaks = avg_random_peaks / true_peak_count if true_peak_count > 0 else 0.0
    fmr_intensity = avg_random_intensity / true_intensity_sum if true_intensity_sum > 0 else 0.0

    return FalseMatchRate(
        fmr_peaks=fmr_peaks,
        fmr_intensity=fmr_intensity,
        matched_peaks=true_peak_count,
        matched_intensity=true_intensity_sum,
        avg_random_peaks=avg_random_peaks,
        avg_random_intensity=avg_random_intensity,
        n_shifts=len(shifts)
    )


def calculate_annotation_statistics(
    matched_ions: List[MatchedIon],
    theoretical_ions: List[TheoreticalIon],
    exp_mz: np.ndarray,
    exp_intensity: np.ndarray,
    peptide_length: int,
    min_sn: float = 5.0
) -> Dict:
    """
    Calculate comprehensive annotation statistics.

    Returns statistics similar to those reported by Annotator:
    - Sequence coverage (backbone bonds fragmented)
    - Peaks annotated (fraction of experimental peaks)
    - Intensity annotated (fraction of total intensity)
    - Fragments found (fraction of theoretical fragments)

    Args:
        matched_ions: List of matched ions
        theoretical_ions: List of all theoretical ions
        exp_mz: Experimental m/z array
        exp_intensity: Experimental intensity array
        peptide_length: Length of the peptide sequence
        min_sn: Minimum signal-to-noise ratio for bond coverage counting.
            Noise is estimated as the median peak intensity.

    Returns:
        Dictionary with annotation statistics
    """
    # Calculate sequence coverage (unique backbone positions fragmented)
    total_bonds = peptide_length - 1
    n_term_positions = set()  # b, c ions
    c_term_positions = set()  # y, z ions

    # Apply S/N filter for bond coverage to prevent noise peaks from
    # inflating coverage statistics
    noise_level = np.median(exp_intensity) if len(exp_intensity) > 0 else 1.0
    min_intensity = noise_level * min_sn

    for ion in matched_ions:
        if ion.exp_intensity < min_intensity:
            continue
        if ion.ion_type in ['b', 'c'] and ion.ion_number > 0:
            n_term_positions.add(ion.ion_number)
        elif ion.ion_type in ['y', 'z'] and ion.ion_number > 0:
            c_term_positions.add(ion.ion_number)

    # A bond is covered if we have either N-term or C-term fragment
    covered_bonds = set()
    for pos in n_term_positions:
        covered_bonds.add(pos)  # b_n/c_n covers bond after residue n
    for pos in c_term_positions:
        covered_bonds.add(peptide_length - pos)  # y_n/z_n covers bond before last n residues

    sequence_coverage = len(covered_bonds) / total_bonds if total_bonds > 0 else 0.0

    # Calculate peaks annotated
    matched_mz_set = set(ion.exp_mz for ion in matched_ions)
    peaks_annotated = len(matched_mz_set) / len(exp_mz) if len(exp_mz) > 0 else 0.0

    # Calculate intensity annotated
    total_intensity = np.sum(exp_intensity)
    matched_intensity = sum(ion.exp_intensity for ion in matched_ions)
    intensity_annotated = matched_intensity / total_intensity if total_intensity > 0 else 0.0

    # Calculate fragments found (theoretical fragments that were matched)
    # Count unique theoretical fragments (by type, number, charge)
    theo_backbone = [(ion.ion_type, ion.ion_number, ion.charge)
                     for ion in theoretical_ions
                     if ion.ion_type in ['b', 'y', 'c', 'z'] and not ion.neutral_loss]
    theo_backbone_unique = set(theo_backbone)

    matched_backbone = [(ion.ion_type, ion.ion_number, ion.charge)
                        for ion in matched_ions
                        if ion.ion_type in ['b', 'y', 'c', 'z'] and not ion.neutral_loss]
    matched_backbone_unique = set(matched_backbone)

    fragments_found = len(matched_backbone_unique) / len(theo_backbone_unique) if len(theo_backbone_unique) > 0 else 0.0

    return {
        'sequence_coverage': sequence_coverage,
        'sequence_coverage_bonds': f"{len(covered_bonds)}/{total_bonds}",
        'peaks_annotated': peaks_annotated,
        'peaks_annotated_count': f"{len(matched_mz_set)}/{len(exp_mz)}",
        'intensity_annotated': intensity_annotated,
        'fragments_found': fragments_found,
        'fragments_found_count': f"{len(matched_backbone_unique)}/{len(theo_backbone_unique)}",
        'n_term_coverage': n_term_positions,
        'c_term_coverage': c_term_positions,
    }


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def parse_modifications_from_string(mod_string: str) -> List[Dict]:
    """
    Parse modification string from FragPipe format.

    Example: "N-term(229.1629),4S(528.2859),19K(229.1629)"
    """
    if not mod_string or mod_string == 'nan':
        return []

    mods = []
    for part in mod_string.split(','):
        part = part.strip()
        if not part or '(' not in part:
            continue

        try:
            mass_str = part[part.find('(')+1:part.find(')')]
            mass = float(mass_str)
            position_part = part[:part.find('(')]

            if position_part == 'N-term':
                mods.append({'position': 0, 'residue': 'N-term', 'mass': mass})
            elif position_part == 'C-term':
                mods.append({'position': -1, 'residue': 'C-term', 'mass': mass})
            else:
                # Parse position and residue
                pos = ''
                res = ''
                for char in position_part:
                    if char.isdigit():
                        pos += char
                    else:
                        res += char
                if pos and res:
                    mods.append({'position': int(pos), 'residue': res, 'mass': mass})
        except:
            continue

    return mods


# =============================================================================
# TEST / DEMO
# =============================================================================

if __name__ == "__main__":
    # Test with the example peptide from your data
    peptide = "AGYSQGATQYTQAQQTR"
    mod_string = "N-term(229.1629),4S(528.2859)"
    precursor_charge = 3

    print("=" * 70)
    print("Fragment Ion Calculator Test")
    print("=" * 70)
    print(f"Peptide: {peptide}")
    print(f"Modifications: {mod_string}")
    print(f"Precursor charge: +{precursor_charge}")

    # Parse modifications
    modifications = parse_modifications_from_string(mod_string)
    print(f"\nParsed modifications:")
    for mod in modifications:
        print(f"  Position {mod['position']}: {mod['residue']} +{mod['mass']:.4f}")

    # Create calculator
    calc = FragmentCalculator(peptide, modifications, precursor_charge)

    print(f"\nPrecursor mass: {calc.precursor_mass:.4f} Da")
    print(f"Precursor m/z: {calc.precursor_mz:.4f} (z={precursor_charge})")

    # Calculate all ions
    all_ions = calc.calculate_all_ions()

    print("\n" + "=" * 70)
    print("Theoretical Ions Summary")
    print("=" * 70)

    for ion_type, ions in all_ions.items():
        if ions:
            print(f"\n{ion_type} ions: {len(ions)}")
            # Show first few
            for ion in ions[:5]:
                print(f"  {ion.annotation}: m/z = {ion.mz:.4f}")
            if len(ions) > 5:
                print(f"  ... and {len(ions) - 5} more")

    # Test with example spectrum data
    print("\n" + "=" * 70)
    print("Peak Matching Test")
    print("=" * 70)

    # Load a test spectrum
    import pandas as pd
    spec_file = "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/point_to_point_response/extracted_spectra_EThcD_HEK293T/Q96KR1_S195_18782.csv"

    try:
        spec_df = pd.read_csv(spec_file)
        exp_mz = spec_df['mz'].values
        exp_intensity = spec_df['intensity'].values

        # Get all theoretical ions
        all_theo_ions = calc.get_all_ions_flat()
        print(f"\nTotal theoretical ions: {len(all_theo_ions)}")

        # Match peaks
        matched = match_peaks(all_theo_ions, exp_mz, exp_intensity, tolerance_ppm=20.0)
        print(f"Matched ions: {len(matched)}")

        # Summarize by ion type
        from collections import Counter
        type_counts = Counter(ion.ion_type for ion in matched)
        print("\nMatched ions by type:")
        for ion_type, count in sorted(type_counts.items()):
            print(f"  {ion_type}: {count}")

        # Show top matched ions by intensity
        print("\nTop 10 matched ions by intensity:")
        matched_sorted = sorted(matched, key=lambda x: x.exp_intensity, reverse=True)
        for i, ion in enumerate(matched_sorted[:10], 1):
            print(f"  {i}. {ion.annotation}: theo={ion.mz:.4f}, exp={ion.exp_mz:.4f}, "
                  f"err={ion.mass_error_ppm:.1f}ppm, int={ion.exp_intensity:.0f}")

        # =================================================================
        # FALSE MATCH RATE CALCULATION
        # =================================================================
        print("\n" + "=" * 70)
        print("False Match Rate Calculation")
        print("(Spectrum shifting method from Schulte et al., Anal. Chem. 2025)")
        print("=" * 70)

        fmr = calculate_false_match_rate(
            all_theo_ions, exp_mz, exp_intensity,
            tolerance_ppm=20.0, shift_range=25.0, shift_step=1.0
        )

        print(f"\nTrue spectrum matches:")
        print(f"  Matched peaks: {fmr.matched_peaks}")
        print(f"  Matched intensity: {fmr.matched_intensity:.0f}")

        print(f"\nShifted spectrum matches (averaged over {fmr.n_shifts} shifts):")
        print(f"  Avg random peaks: {fmr.avg_random_peaks:.2f}")
        print(f"  Avg random intensity: {fmr.avg_random_intensity:.0f}")

        print(f"\nFalse Match Rates:")
        print(f"  FMR (peaks): {fmr.fmr_peaks*100:.1f}%")
        print(f"  FMR (intensity): {fmr.fmr_intensity*100:.1f}%")

        print("\nInterpretation:")
        print(f"  ~{fmr.fmr_peaks*100:.1f}% of matched peaks may be spurious")
        print(f"  ~{fmr.fmr_intensity*100:.1f}% of matched intensity may be spurious")

        # =================================================================
        # ANNOTATION STATISTICS
        # =================================================================
        print("\n" + "=" * 70)
        print("Annotation Statistics")
        print("=" * 70)

        stats = calculate_annotation_statistics(
            matched, all_theo_ions, exp_mz, exp_intensity, len(peptide)
        )

        print(f"\nSequence coverage: {stats['sequence_coverage_bonds']} "
              f"({stats['sequence_coverage']*100:.1f}%)")
        print(f"Peaks annotated: {stats['peaks_annotated_count']} "
              f"({stats['peaks_annotated']*100:.1f}%)")
        print(f"Intensity annotated: {stats['intensity_annotated']*100:.1f}%")
        print(f"Fragments found: {stats['fragments_found_count']} "
              f"({stats['fragments_found']*100:.1f}%)")

    except FileNotFoundError:
        print("Test spectrum file not found. Skipping peak matching test.")

    print("\n" + "=" * 70)
    print("Test complete!")
    print("=" * 70)
