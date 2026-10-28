# ABOUTME: The Mechanisms module: reactions and mechanisms as values (raw input, canonical
# ABOUTME: graph, Sig), their naming rules, and the DSL's unchecked reader and printer.
module Mechanisms

import Unicode
using SHA: sha256

const API = (:AbstractEnzymeMechanism, :AllostericEnzymeMechanism, :AllostericMechanism,
    :AllostericRegulator, :CompetitiveInhibitor, :EnzymeMechanism, :EnzymeReaction,
    :Mechanism, :Met, :Metabolites, :RawAllostericMechanism, :RawForm, :RawMechanism,
    :RawSite, :RawStep, :RegulatorySite, :_AnyMech, :_Column, :_canonicalize,
    :_check_aligned, :_check_uniquely_decodable, :_chemistry_groups, :_concrete,
    :_edge_constants, :_excluded_by, :_excludes, :_exclusion_components, :_flat_edges,
    :_form_name, :_form_text, :_group_members, :_groups_with, :_has_inactive_columns,
    :_has_unbalanced_cycle, :_in_active_law, :_is_iso, :_kept_members, :_kreg_params,
    :_metabolite_names, :_number_groups, :_ordered_columns, :_parse_reaction_block,
    :_raw_mechanism, :_reaction_expr, :_reg_K, :_runs_toward_products, :_segment_walk,
    :_sides, :_solve_tags, :_step_text, :allo_states, :allowed_catalytic_multiplicities,
    :atoms, :bound, :cat_allo_states, :catalytic, :catalytic_multiplicity,
    :compile_mechanism, :conformation, :conformation_steps, :copies, :ends, :given_off,
    :group, :has_turnover_cycle, :is_binding, :is_rapid_equilibrium, :ligands,
    :mechanism_id, :members, :metabolite, :metabolite_names, :metabolites, :modifier_copies,
    :modifiers, :multiplicity, :name, :ngroups, :present_species, :products,
    :rapid_equilibrium_steps, :re_segments, :reactants, :reg_type, :regroup, :regulator,
    :regulators, :regulatory_sites, :residual, :segment_graph, :shared_catalytic_site,
    :species, :step_kind, :steps, :substrates, :tag, :taken_up, :with_flags, :with_tags)

include("types.jl")
include("accessors.jl")
include("write.jl")
include("names.jl")
include("dsl.jl")

end
