# ABOUTME: The Mechanisms module: reactions and mechanisms as values (a canonical graph
# ABOUTME: built from nodes and steps, Sig), their naming rules, and their DSL text.
module Mechanisms

import Unicode

const API = (:AbstractEnzymeMechanism, :AllostericEnzymeMechanism, :AllostericMechanism,
    :EnzymeMechanism, :EnzymeReaction, :Mechanism, :Met, :Node, :Reactant, :Regulator,
    :RegulatorySite, :Step, :catalytic, :check_conformation_label, :check_decodable_names,
    :compile_mechanism, :form_text, :groups, :metabolites, :modifiers,
    :runs_toward_products, :step_kind, :step_sides, :step_text)

include("types.jl")
include("accessors.jl")
include("names.jl")
include("text.jl")

end
