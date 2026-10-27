# ABOUTME: The Mechanisms module: reactions and mechanisms as values (raw input, canonical
# ABOUTME: graph, Sig), their naming rules, and the DSL's unchecked reader and printer.
module Mechanisms

import Unicode
using SHA: sha256

const API = (:AbstractEnzymeMechanism, :AllostericEnzymeMechanism, :AllostericMechanism,
    :AllostericRegulator, :CompetitiveInhibitor, :Counts, :Edge, :EnzymeMechanism,
    :EnzymeReaction, :Form, :Mechanism, :Met, :RawAllostericMechanism, :RawForm,
    :RawMechanism, :RawStep, :RegulatorySite, :_AnyMech, :_Column, :_active_law,
    :_active_mask, :_bound_entry, :_build, :_canonicalize, :_check_aligned,
    :_check_uniquely_decodable, :_chemistry_groups, :_concrete, :_edge_constants, :_entry,
    :_entry_names, :_excluded_by, :_excludes, :_exclusion_components, :_flags, :_flat_edges,
    :_form_name, :_form_text, :_group_members, :_groups_with, :_has_inactive_columns,
    :_has_residual, :_has_unbalanced_cycle, :_in_active_law, :_inactive_mask, :_is_binding,
    :_is_chemistry, :_is_iso, :_is_plain_binding, :_is_reactant, :_kept_members,
    :_kreg_params, :_law_mask, :_metabolite_names, :_modifier_mask, :_name_index,
    :_names_mask, :_number_groups, :_ordered_columns, :_parse_reaction_block, :_raw,
    :_raw_form, :_raw_mechanism, :_raw_sites, :_reaction_expr, :_reg_K, :_role_masks,
    :_segments, :_side_score, :_sides, :_solve_tags, :_step_text, :_with, :allo_states,
    :allowed_catalytic_multiplicities, :atoms, :cat_allo_state, :cat_allo_states,
    :catalytic, :catalytic_multiplicity, :compile_mechanism, :copies, :entries, :ligands,
    :mechanism_id, :metabolite, :metabolites, :multiplicity, :name, :ngroups, :products,
    :reactants, :reg_type, :regulator, :regulators, :regulatory_sites,
    :shared_catalytic_site, :substrates)

include("types.jl")
include("names.jl")
include("dsl.jl")

end
