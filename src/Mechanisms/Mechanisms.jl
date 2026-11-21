# ABOUTME: The Mechanisms module: reactions and mechanisms as values (a catalytic scheme in
# ABOUTME: canonical form and its allostery) with the facts they store, the segment and
# ABOUTME: cycle analyses, their naming rules, and their DSL text.
"""Stores reactions and mechanisms as plain data."""
module Mechanisms

import Unicode

public Allostery, CatalyticScheme, EnzymeReaction, Mechanism, Metabolite, EnzymeForm,
    Reactant, Regulator, RegulatorySite, Step, check_conformation_label, form_text,
    has_turnover_cycle, has_unbalanced_cycle, metabolites, rapid_equilibrium_segments,
    step_text, steps_in_active_conformation, steps_in_inactive_conformation

include("types.jl")
include("accessors.jl")
include("names.jl")
include("text.jl")
include("segments.jl")

end
