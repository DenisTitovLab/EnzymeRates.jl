# ABOUTME: The Mechanisms module: reactions and mechanisms as values (a canonical graph
# ABOUTME: built from nodes and steps, Sig), their naming rules, and their DSL text.
"""Stores reactions and mechanisms as plain data."""
module Mechanisms

import Unicode

const API = (:AbstractEnzymeMechanism, :AllostericEnzymeMechanism, :AllostericMechanism,
    :EnzymeMechanism, :EnzymeReaction, :Mechanism, :Metabolite, :EnzymeForm, :Reactant,
    :Regulator, :RegulatorySite, :Step, :catalytic, :check_conformation_label,
    :check_decodable_names, :compile_mechanism, :excludes, :form_text, :has_turnover_cycle,
    :has_unbalanced_cycle, :metabolites, :present_species, :rapid_equilibrium_segments,
    :segment_walk, :step_text)

include("types.jl")
include("accessors.jl")
include("names.jl")
include("text.jl")
include("segments.jl")

end
