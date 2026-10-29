# ABOUTME: The Mechanisms module: reactions and mechanisms as values (raw input, canonical
# ABOUTME: graph, Sig), their naming rules, and the printer of their DSL text.
module Mechanisms

import Unicode

const API = (:AbstractEnzymeMechanism, :AllostericEnzymeMechanism, :AllostericMechanism,
    :EnzymeMechanism, :EnzymeReaction, :Mechanism, :Met, :Metabolites, :Node,
    :RawAllostericMechanism, :RawForm, :RawMechanism, :RawSite, :RawStep, :Reactant,
    :Regulator, :RegulatorySite, :bound, :canonical_graph,
    :catalytic, :catalytic_multiplicity, :check_conformation_label,
    :check_decodable_names, :compile_mechanism, :conformation, :conformation_steps,
    :constant_columns, :copies, :ends, :excludes, :exclusion_components, :form_text,
    :given_off, :group, :groups, :has_turnover_cycle, :has_unbalanced_cycle,
    :is_rapid_equilibrium, :label_collision_error, :ligands, :members, :metabolite_names,
    :metabolites, :modifier_copies, :modifiers, :multiplicity, :mwc_tags, :name, :ngroups,
    :present_species, :products, :rapid_equilibrium_steps, :re_segments, :regroup,
    :regulator_constant, :regulatory_sites, :residual, :runs_toward_products,
    :segment_graph, :segment_walk, :species, :Step, :step_kind, :step_text, :steps,
    :substrates, :tag, :taken_up, :with_flags, :with_tags)

include("types.jl")
include("accessors.jl")
include("write.jl")
include("names.jl")
include("text.jl")

end
