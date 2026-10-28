# ABOUTME: The Mechanisms module: reactions and mechanisms as values (raw input, canonical
# ABOUTME: graph, Sig), their naming rules, and the DSL's unchecked reader and printer.
module Mechanisms

import Unicode

const API = (:AbstractEnzymeMechanism, :AllostericEnzymeMechanism, :AllostericMechanism,
    :AllostericRegulator, :CompetitiveInhibitor, :EnzymeMechanism, :EnzymeReaction,
    :Mechanism, :Met, :Metabolites, :RawAllostericMechanism, :RawForm, :RawMechanism,
    :RawSite, :RawStep, :RegulatorySite, :allo_states, :allowed_catalytic_multiplicities,
    :atoms, :bound, :canonical_graph, :cat_allo_states, :catalytic, :catalytic_multiplicity,
    :check_uniquely_decodable, :compile_mechanism, :conformation, :conformation_steps,
    :constant_columns, :copies, :ends, :excludes, :exclusion_components, :form_text,
    :given_off, :group, :has_turnover_cycle, :has_unbalanced_cycle, :is_binding,
    :is_rapid_equilibrium, :ligands, :members, :metabolite,
    :metabolite_names, :metabolites, :modifier_copies, :modifiers, :multiplicity, :name,
    :ngroups, :present_species, :products, :rapid_equilibrium_steps, :raw_mechanism,
    :re_segments, :reactants, :reaction_expr, :reg_type, :regroup, :regulator,
    :regulator_constant, :regulators, :regulatory_sites, :residual, :runs_toward_products,
    :segment_graph, :segment_walk, :shared_catalytic_site, :species, :step_kind, :step_text,
    :steps, :substrates, :tag, :taken_up, :with_flags, :with_tags)

include("types.jl")
include("accessors.jl")
include("write.jl")
include("names.jl")
include("dsl.jl")

end
