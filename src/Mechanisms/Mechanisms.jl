# ABOUTME: The Mechanisms module: reactions and mechanisms as values (raw input, canonical
# ABOUTME: graph, Sig), their naming rules, and the printer of their DSL text.
module Mechanisms

import Unicode

const API = (:AbstractEnzymeMechanism, :AllostericEnzymeMechanism, :AllostericMechanism,
    :EnzymeMechanism, :EnzymeReaction, :Mechanism, :Met, :Node, :RawAllostericMechanism,
    :RawForm, :RawMechanism, :RawSite, :RawStep, :Reactant, :Regulator, :RegulatorySite,
    :Step, :canonical_graph, :catalytic, :check_conformation_label, :check_decodable_names,
    :compile_mechanism, :form_text, :groups, :label_collision_error, :metabolites,
    :modifiers, :regroup, :runs_toward_products, :step_kind, :step_sides, :step_text,
    :with_flags, :with_tags)

include("types.jl")
include("accessors.jl")
include("write.jl")
include("names.jl")
include("text.jl")

end
