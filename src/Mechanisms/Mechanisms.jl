# ABOUTME: The Mechanisms module: reactions and mechanisms as values (raw input, canonical
# ABOUTME: graph, Sig), their naming rules, and the printer of their DSL text.
module Mechanisms

import Unicode

const API = (:AbstractEnzymeMechanism, :AllostericEnzymeMechanism, :AllostericMechanism,
    :EnzymeMechanism, :EnzymeReaction, :Mechanism, :Met, :Metabolites,
    :RawAllostericMechanism, :RawForm, :RawMechanism, :RawSite, :RawStep, :Reactant,
    :Regulator, :RegulatorySite, :allo_states, :bound, :canonical_graph, :cat_allo_states,
    :catalytic, :catalytic_multiplicity, :check_conformation_label,
    :check_decodable_names, :compile_mechanism, :conformation, :conformation_steps,
    :constant_columns, :copies, :ends, :excludes, :exclusion_components, :form_text,
    :given_off, :group, :has_turnover_cycle, :has_unbalanced_cycle, :is_binding,
    :is_rapid_equilibrium, :label_collision_error, :ligands, :members, :metabolite_names,
    :metabolites, :modifier_copies, :modifiers, :multiplicity, :name, :ngroups,
    :present_species, :products, :rapid_equilibrium_steps, :re_segments, :regroup,
    :regulator_constant, :regulatory_sites, :residual, :runs_toward_products,
    :segment_graph, :segment_walk, :species, :step_kind, :step_text, :steps, :substrates,
    :tag, :taken_up, :with_flags, :with_tags)

include("types.jl")
include("accessors.jl")
include("write.jl")
include("names.jl")
include("text.jl")

end
