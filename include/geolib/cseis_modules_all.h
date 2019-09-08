#ifndef CS_MODULES_ALL_H
#define CS_MODULES_ALL_H

namespace cseis_system {
  class csTrace;
  class csParamManager;
  class csTraceGather;
  class csParamDef;
  class csInitPhaseEnv;
  class csExecPhaseEnv;
  class csLogWriter;
}

const int N_METHODS_SINGLE = 47;
const int N_METHODS_MULTI  = 36;
const int N_METHODS = N_METHODS_MULTI + N_METHODS_SINGLE;

void init_mod_attribute_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_beam_forming_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_bin_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_ccp_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_cmp_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_concatenate_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_convolution_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_correlation_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_debias_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_designature_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_despike_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_else_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_elseif_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_endif_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_endsplit_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_ens_define_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_fft_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_fft_2d_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_filter_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_fxdecon_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_gain_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_geotools_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_hdr_del_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_hdr_math_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_hdr_math_ens_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_hdr_print_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_hdr_set_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_histogram_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_hodogram_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_if_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_image_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_input_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_input_ascii_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_input_create_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_input_rsf_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_input_segd_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_input_segy_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_input_sinewave_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_kill_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_kill_ens_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_lmo_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_mirror_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_mute_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_nmo_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_off2angle_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_orient_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_orient_convert_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_output_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_output_rsf_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_output_segy_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_overlap_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_p190_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_picking_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_poscalc_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_pz_sum_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_ray2d_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_read_ascii_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_repeat_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_resample_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_resequence_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_rms_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_rotate_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_scaling_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_select_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_select_time_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_semblance_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_sort_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_split_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_splitting_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_stack_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_statics_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_sumodule_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_test_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_test_multi_ensemble_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_test_multi_fixed_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_time_slice_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_time_stretch_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_trc_add_ens_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_trc_interpol_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_trc_math_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_trc_math_ens_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_trc_print_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);
void init_mod_trc_split_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);

bool exec_mod_dummy_single_(cseis_system::csTrace*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*) { return true; }
bool exec_mod_attribute_(cseis_system::csTrace*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
bool exec_mod_bin_(cseis_system::csTrace*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
bool exec_mod_cmp_(cseis_system::csTrace*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
bool exec_mod_convolution_(cseis_system::csTrace*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
bool exec_mod_designature_(cseis_system::csTrace*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
bool exec_mod_despike_(cseis_system::csTrace*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
bool exec_mod_else_(cseis_system::csTrace*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
bool exec_mod_elseif_(cseis_system::csTrace*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
bool exec_mod_endif_(cseis_system::csTrace*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
bool exec_mod_endsplit_(cseis_system::csTrace*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
bool exec_mod_ens_define_(cseis_system::csTrace*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
bool exec_mod_fft_(cseis_system::csTrace*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
bool exec_mod_gain_(cseis_system::csTrace*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
bool exec_mod_geotools_(cseis_system::csTrace*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
bool exec_mod_hdr_del_(cseis_system::csTrace*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
bool exec_mod_hdr_math_(cseis_system::csTrace*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
bool exec_mod_hdr_print_(cseis_system::csTrace*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
bool exec_mod_hdr_set_(cseis_system::csTrace*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
bool exec_mod_if_(cseis_system::csTrace*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
bool exec_mod_image_(cseis_system::csTrace*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
bool exec_mod_input_(cseis_system::csTrace*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
bool exec_mod_input_ascii_(cseis_system::csTrace*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
bool exec_mod_input_create_(cseis_system::csTrace*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
bool exec_mod_input_rsf_(cseis_system::csTrace*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
bool exec_mod_input_segd_(cseis_system::csTrace*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
bool exec_mod_input_segy_(cseis_system::csTrace*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
bool exec_mod_input_sinewave_(cseis_system::csTrace*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
bool exec_mod_kill_(cseis_system::csTrace*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
bool exec_mod_lmo_(cseis_system::csTrace*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
bool exec_mod_mirror_(cseis_system::csTrace*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
bool exec_mod_mute_(cseis_system::csTrace*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
bool exec_mod_nmo_(cseis_system::csTrace*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
bool exec_mod_output_(cseis_system::csTrace*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
bool exec_mod_output_rsf_(cseis_system::csTrace*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
bool exec_mod_output_segy_(cseis_system::csTrace*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
bool exec_mod_p190_(cseis_system::csTrace*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
bool exec_mod_read_ascii_(cseis_system::csTrace*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
bool exec_mod_resample_(cseis_system::csTrace*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
bool exec_mod_rms_(cseis_system::csTrace*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
bool exec_mod_scaling_(cseis_system::csTrace*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
bool exec_mod_select_(cseis_system::csTrace*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
bool exec_mod_select_time_(cseis_system::csTrace*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
bool exec_mod_statics_(cseis_system::csTrace*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
bool exec_mod_test_(cseis_system::csTrace*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
bool exec_mod_time_stretch_(cseis_system::csTrace*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
bool exec_mod_trc_math_(cseis_system::csTrace*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
bool exec_mod_trc_print_(cseis_system::csTrace*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);

void exec_mod_dummy_multi_(cseis_system::csTraceGather*, int*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*) {}
void exec_mod_beam_forming_(cseis_system::csTraceGather*, int*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
void exec_mod_ccp_(cseis_system::csTraceGather*, int*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
void exec_mod_concatenate_(cseis_system::csTraceGather*, int*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
void exec_mod_correlation_(cseis_system::csTraceGather*, int*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
void exec_mod_debias_(cseis_system::csTraceGather*, int*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
void exec_mod_fft_2d_(cseis_system::csTraceGather*, int*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
void exec_mod_filter_(cseis_system::csTraceGather*, int*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
void exec_mod_fxdecon_(cseis_system::csTraceGather*, int*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
void exec_mod_hdr_math_ens_(cseis_system::csTraceGather*, int*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
void exec_mod_histogram_(cseis_system::csTraceGather*, int*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
void exec_mod_hodogram_(cseis_system::csTraceGather*, int*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
void exec_mod_kill_ens_(cseis_system::csTraceGather*, int*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
void exec_mod_off2angle_(cseis_system::csTraceGather*, int*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
void exec_mod_orient_(cseis_system::csTraceGather*, int*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
void exec_mod_orient_convert_(cseis_system::csTraceGather*, int*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
void exec_mod_overlap_(cseis_system::csTraceGather*, int*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
void exec_mod_picking_(cseis_system::csTraceGather*, int*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
void exec_mod_poscalc_(cseis_system::csTraceGather*, int*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
void exec_mod_pz_sum_(cseis_system::csTraceGather*, int*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
void exec_mod_ray2d_(cseis_system::csTraceGather*, int*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
void exec_mod_repeat_(cseis_system::csTraceGather*, int*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
void exec_mod_resequence_(cseis_system::csTraceGather*, int*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
void exec_mod_rotate_(cseis_system::csTraceGather*, int*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
void exec_mod_semblance_(cseis_system::csTraceGather*, int*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
void exec_mod_sort_(cseis_system::csTraceGather*, int*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
void exec_mod_split_(cseis_system::csTraceGather*, int*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
void exec_mod_splitting_(cseis_system::csTraceGather*, int*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
void exec_mod_stack_(cseis_system::csTraceGather*, int*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
void exec_mod_sumodule_(cseis_system::csTraceGather*, int*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
void exec_mod_test_multi_ensemble_(cseis_system::csTraceGather*, int*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
void exec_mod_test_multi_fixed_(cseis_system::csTraceGather*, int*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
void exec_mod_time_slice_(cseis_system::csTraceGather*, int*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
void exec_mod_trc_add_ens_(cseis_system::csTraceGather*, int*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
void exec_mod_trc_interpol_(cseis_system::csTraceGather*, int*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
void exec_mod_trc_math_ens_(cseis_system::csTraceGather*, int*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);
void exec_mod_trc_split_(cseis_system::csTraceGather*, int*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);

void(*METHODS_INIT[N_METHODS])( cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter* ) = {
  init_mod_attribute_,
  init_mod_bin_,
  init_mod_cmp_,
  init_mod_convolution_,
  init_mod_designature_,
  init_mod_despike_,
  init_mod_else_,
  init_mod_elseif_,
  init_mod_endif_,
  init_mod_endsplit_,
  init_mod_ens_define_,
  init_mod_fft_,
  init_mod_gain_,
  init_mod_geotools_,
  init_mod_hdr_del_,
  init_mod_hdr_math_,
  init_mod_hdr_print_,
  init_mod_hdr_set_,
  init_mod_if_,
  init_mod_image_,
  init_mod_input_,
  init_mod_input_ascii_,
  init_mod_input_create_,
  init_mod_input_rsf_,
  init_mod_input_segd_,
  init_mod_input_segy_,
  init_mod_input_sinewave_,
  init_mod_kill_,
  init_mod_lmo_,
  init_mod_mirror_,
  init_mod_mute_,
  init_mod_nmo_,
  init_mod_output_,
  init_mod_output_rsf_,
  init_mod_output_segy_,
  init_mod_p190_,
  init_mod_read_ascii_,
  init_mod_resample_,
  init_mod_rms_,
  init_mod_scaling_,
  init_mod_select_,
  init_mod_select_time_,
  init_mod_statics_,
  init_mod_test_,
  init_mod_time_stretch_,
  init_mod_trc_math_,
  init_mod_trc_print_,
  init_mod_beam_forming_,
  init_mod_ccp_,
  init_mod_concatenate_,
  init_mod_correlation_,
  init_mod_debias_,
  init_mod_fft_2d_,
  init_mod_filter_,
  init_mod_fxdecon_,
  init_mod_hdr_math_ens_,
  init_mod_histogram_,
  init_mod_hodogram_,
  init_mod_kill_ens_,
  init_mod_off2angle_,
  init_mod_orient_,
  init_mod_orient_convert_,
  init_mod_overlap_,
  init_mod_picking_,
  init_mod_poscalc_,
  init_mod_pz_sum_,
  init_mod_ray2d_,
  init_mod_repeat_,
  init_mod_resequence_,
  init_mod_rotate_,
  init_mod_semblance_,
  init_mod_sort_,
  init_mod_split_,
  init_mod_splitting_,
  init_mod_stack_,
  init_mod_sumodule_,
  init_mod_test_multi_ensemble_,
  init_mod_test_multi_fixed_,
  init_mod_time_slice_,
  init_mod_trc_add_ens_,
  init_mod_trc_interpol_,
  init_mod_trc_math_ens_,
  init_mod_trc_split_
};

bool(*METHODS_EXEC_SINGLE[N_METHODS])(cseis_system::csTrace*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*) = {
  exec_mod_attribute_,
  exec_mod_bin_,
  exec_mod_cmp_,
  exec_mod_convolution_,
  exec_mod_designature_,
  exec_mod_despike_,
  exec_mod_else_,
  exec_mod_elseif_,
  exec_mod_endif_,
  exec_mod_endsplit_,
  exec_mod_ens_define_,
  exec_mod_fft_,
  exec_mod_gain_,
  exec_mod_geotools_,
  exec_mod_hdr_del_,
  exec_mod_hdr_math_,
  exec_mod_hdr_print_,
  exec_mod_hdr_set_,
  exec_mod_if_,
  exec_mod_image_,
  exec_mod_input_,
  exec_mod_input_ascii_,
  exec_mod_input_create_,
  exec_mod_input_rsf_,
  exec_mod_input_segd_,
  exec_mod_input_segy_,
  exec_mod_input_sinewave_,
  exec_mod_kill_,
  exec_mod_lmo_,
  exec_mod_mirror_,
  exec_mod_mute_,
  exec_mod_nmo_,
  exec_mod_output_,
  exec_mod_output_rsf_,
  exec_mod_output_segy_,
  exec_mod_p190_,
  exec_mod_read_ascii_,
  exec_mod_resample_,
  exec_mod_rms_,
  exec_mod_scaling_,
  exec_mod_select_,
  exec_mod_select_time_,
  exec_mod_statics_,
  exec_mod_test_,
  exec_mod_time_stretch_,
  exec_mod_trc_math_,
  exec_mod_trc_print_,
  exec_mod_dummy_single_,
  exec_mod_dummy_single_,
  exec_mod_dummy_single_,
  exec_mod_dummy_single_,
  exec_mod_dummy_single_,
  exec_mod_dummy_single_,
  exec_mod_dummy_single_,
  exec_mod_dummy_single_,
  exec_mod_dummy_single_,
  exec_mod_dummy_single_,
  exec_mod_dummy_single_,
  exec_mod_dummy_single_,
  exec_mod_dummy_single_,
  exec_mod_dummy_single_,
  exec_mod_dummy_single_,
  exec_mod_dummy_single_,
  exec_mod_dummy_single_,
  exec_mod_dummy_single_,
  exec_mod_dummy_single_,
  exec_mod_dummy_single_,
  exec_mod_dummy_single_,
  exec_mod_dummy_single_,
  exec_mod_dummy_single_,
  exec_mod_dummy_single_,
  exec_mod_dummy_single_,
  exec_mod_dummy_single_,
  exec_mod_dummy_single_,
  exec_mod_dummy_single_,
  exec_mod_dummy_single_,
  exec_mod_dummy_single_,
  exec_mod_dummy_single_,
  exec_mod_dummy_single_,
  exec_mod_dummy_single_,
  exec_mod_dummy_single_,
  exec_mod_dummy_single_,
  exec_mod_dummy_single_
};

void(*METHODS_EXEC_MULTI[N_METHODS])(cseis_system::csTraceGather*, int*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*) = {
  exec_mod_dummy_multi_,
  exec_mod_dummy_multi_,
  exec_mod_dummy_multi_,
  exec_mod_dummy_multi_,
  exec_mod_dummy_multi_,
  exec_mod_dummy_multi_,
  exec_mod_dummy_multi_,
  exec_mod_dummy_multi_,
  exec_mod_dummy_multi_,
  exec_mod_dummy_multi_,
  exec_mod_dummy_multi_,
  exec_mod_dummy_multi_,
  exec_mod_dummy_multi_,
  exec_mod_dummy_multi_,
  exec_mod_dummy_multi_,
  exec_mod_dummy_multi_,
  exec_mod_dummy_multi_,
  exec_mod_dummy_multi_,
  exec_mod_dummy_multi_,
  exec_mod_dummy_multi_,
  exec_mod_dummy_multi_,
  exec_mod_dummy_multi_,
  exec_mod_dummy_multi_,
  exec_mod_dummy_multi_,
  exec_mod_dummy_multi_,
  exec_mod_dummy_multi_,
  exec_mod_dummy_multi_,
  exec_mod_dummy_multi_,
  exec_mod_dummy_multi_,
  exec_mod_dummy_multi_,
  exec_mod_dummy_multi_,
  exec_mod_dummy_multi_,
  exec_mod_dummy_multi_,
  exec_mod_dummy_multi_,
  exec_mod_dummy_multi_,
  exec_mod_dummy_multi_,
  exec_mod_dummy_multi_,
  exec_mod_dummy_multi_,
  exec_mod_dummy_multi_,
  exec_mod_dummy_multi_,
  exec_mod_dummy_multi_,
  exec_mod_dummy_multi_,
  exec_mod_dummy_multi_,
  exec_mod_dummy_multi_,
  exec_mod_dummy_multi_,
  exec_mod_dummy_multi_,
  exec_mod_dummy_multi_,
  exec_mod_beam_forming_,
  exec_mod_ccp_,
  exec_mod_concatenate_,
  exec_mod_correlation_,
  exec_mod_debias_,
  exec_mod_fft_2d_,
  exec_mod_filter_,
  exec_mod_fxdecon_,
  exec_mod_hdr_math_ens_,
  exec_mod_histogram_,
  exec_mod_hodogram_,
  exec_mod_kill_ens_,
  exec_mod_off2angle_,
  exec_mod_orient_,
  exec_mod_orient_convert_,
  exec_mod_overlap_,
  exec_mod_picking_,
  exec_mod_poscalc_,
  exec_mod_pz_sum_,
  exec_mod_ray2d_,
  exec_mod_repeat_,
  exec_mod_resequence_,
  exec_mod_rotate_,
  exec_mod_semblance_,
  exec_mod_sort_,
  exec_mod_split_,
  exec_mod_splitting_,
  exec_mod_stack_,
  exec_mod_sumodule_,
  exec_mod_test_multi_ensemble_,
  exec_mod_test_multi_fixed_,
  exec_mod_time_slice_,
  exec_mod_trc_add_ens_,
  exec_mod_trc_interpol_,
  exec_mod_trc_math_ens_,
  exec_mod_trc_split_
};

std::string NAMES[N_METHODS] = {
  std::string("ATTRIBUTE"),
  std::string("BIN"),
  std::string("CMP"),
  std::string("CONVOLUTION"),
  std::string("DESIGNATURE"),
  std::string("DESPIKE"),
  std::string("ELSE"),
  std::string("ELSEIF"),
  std::string("ENDIF"),
  std::string("ENDSPLIT"),
  std::string("ENS_DEFINE"),
  std::string("FFT"),
  std::string("GAIN"),
  std::string("GEOTOOLS"),
  std::string("HDR_DEL"),
  std::string("HDR_MATH"),
  std::string("HDR_PRINT"),
  std::string("HDR_SET"),
  std::string("IF"),
  std::string("IMAGE"),
  std::string("INPUT"),
  std::string("INPUT_ASCII"),
  std::string("INPUT_CREATE"),
  std::string("INPUT_RSF"),
  std::string("INPUT_SEGD"),
  std::string("INPUT_SEGY"),
  std::string("INPUT_SINEWAVE"),
  std::string("KILL"),
  std::string("LMO"),
  std::string("MIRROR"),
  std::string("MUTE"),
  std::string("NMO"),
  std::string("OUTPUT"),
  std::string("OUTPUT_RSF"),
  std::string("OUTPUT_SEGY"),
  std::string("P190"),
  std::string("READ_ASCII"),
  std::string("RESAMPLE"),
  std::string("RMS"),
  std::string("SCALING"),
  std::string("SELECT"),
  std::string("SELECT_TIME"),
  std::string("STATICS"),
  std::string("TEST"),
  std::string("TIME_STRETCH"),
  std::string("TRC_MATH"),
  std::string("TRC_PRINT"),
  std::string("BEAM_FORMING"),
  std::string("CCP"),
  std::string("CONCATENATE"),
  std::string("CORRELATION"),
  std::string("DEBIAS"),
  std::string("FFT_2D"),
  std::string("FILTER"),
  std::string("FXDECON"),
  std::string("HDR_MATH_ENS"),
  std::string("HISTOGRAM"),
  std::string("HODOGRAM"),
  std::string("KILL_ENS"),
  std::string("OFF2ANGLE"),
  std::string("ORIENT"),
  std::string("ORIENT_CONVERT"),
  std::string("OVERLAP"),
  std::string("PICKING"),
  std::string("POSCALC"),
  std::string("PZ_SUM"),
  std::string("RAY2D"),
  std::string("REPEAT"),
  std::string("RESEQUENCE"),
  std::string("ROTATE"),
  std::string("SEMBLANCE"),
  std::string("SORT"),
  std::string("SPLIT"),
  std::string("SPLITTING"),
  std::string("STACK"),
  std::string("SUMODULE"),
  std::string("TEST_MULTI_ENSEMBLE"),
  std::string("TEST_MULTI_FIXED"),
  std::string("TIME_SLICE"),
  std::string("TRC_ADD_ENS"),
  std::string("TRC_INTERPOL"),
  std::string("TRC_MATH_ENS"),
  std::string("TRC_SPLIT")
};

void params_mod_attribute_(cseis_system::csParamDef*);
void params_mod_beam_forming_(cseis_system::csParamDef*);
void params_mod_bin_(cseis_system::csParamDef*);
void params_mod_ccp_(cseis_system::csParamDef*);
void params_mod_cmp_(cseis_system::csParamDef*);
void params_mod_concatenate_(cseis_system::csParamDef*);
void params_mod_convolution_(cseis_system::csParamDef*);
void params_mod_correlation_(cseis_system::csParamDef*);
void params_mod_debias_(cseis_system::csParamDef*);
void params_mod_designature_(cseis_system::csParamDef*);
void params_mod_despike_(cseis_system::csParamDef*);
void params_mod_else_(cseis_system::csParamDef*);
void params_mod_elseif_(cseis_system::csParamDef*);
void params_mod_endif_(cseis_system::csParamDef*);
void params_mod_endsplit_(cseis_system::csParamDef*);
void params_mod_ens_define_(cseis_system::csParamDef*);
void params_mod_fft_(cseis_system::csParamDef*);
void params_mod_fft_2d_(cseis_system::csParamDef*);
void params_mod_filter_(cseis_system::csParamDef*);
void params_mod_fxdecon_(cseis_system::csParamDef*);
void params_mod_gain_(cseis_system::csParamDef*);
void params_mod_geotools_(cseis_system::csParamDef*);
void params_mod_hdr_del_(cseis_system::csParamDef*);
void params_mod_hdr_math_(cseis_system::csParamDef*);
void params_mod_hdr_math_ens_(cseis_system::csParamDef*);
void params_mod_hdr_print_(cseis_system::csParamDef*);
void params_mod_hdr_set_(cseis_system::csParamDef*);
void params_mod_histogram_(cseis_system::csParamDef*);
void params_mod_hodogram_(cseis_system::csParamDef*);
void params_mod_if_(cseis_system::csParamDef*);
void params_mod_image_(cseis_system::csParamDef*);
void params_mod_input_(cseis_system::csParamDef*);
void params_mod_input_ascii_(cseis_system::csParamDef*);
void params_mod_input_create_(cseis_system::csParamDef*);
void params_mod_input_rsf_(cseis_system::csParamDef*);
void params_mod_input_segd_(cseis_system::csParamDef*);
void params_mod_input_segy_(cseis_system::csParamDef*);
void params_mod_input_sinewave_(cseis_system::csParamDef*);
void params_mod_kill_(cseis_system::csParamDef*);
void params_mod_kill_ens_(cseis_system::csParamDef*);
void params_mod_lmo_(cseis_system::csParamDef*);
void params_mod_mirror_(cseis_system::csParamDef*);
void params_mod_mute_(cseis_system::csParamDef*);
void params_mod_nmo_(cseis_system::csParamDef*);
void params_mod_off2angle_(cseis_system::csParamDef*);
void params_mod_orient_(cseis_system::csParamDef*);
void params_mod_orient_convert_(cseis_system::csParamDef*);
void params_mod_output_(cseis_system::csParamDef*);
void params_mod_output_rsf_(cseis_system::csParamDef*);
void params_mod_output_segy_(cseis_system::csParamDef*);
void params_mod_overlap_(cseis_system::csParamDef*);
void params_mod_p190_(cseis_system::csParamDef*);
void params_mod_picking_(cseis_system::csParamDef*);
void params_mod_poscalc_(cseis_system::csParamDef*);
void params_mod_pz_sum_(cseis_system::csParamDef*);
void params_mod_ray2d_(cseis_system::csParamDef*);
void params_mod_read_ascii_(cseis_system::csParamDef*);
void params_mod_repeat_(cseis_system::csParamDef*);
void params_mod_resample_(cseis_system::csParamDef*);
void params_mod_resequence_(cseis_system::csParamDef*);
void params_mod_rms_(cseis_system::csParamDef*);
void params_mod_rotate_(cseis_system::csParamDef*);
void params_mod_scaling_(cseis_system::csParamDef*);
void params_mod_select_(cseis_system::csParamDef*);
void params_mod_select_time_(cseis_system::csParamDef*);
void params_mod_semblance_(cseis_system::csParamDef*);
void params_mod_sort_(cseis_system::csParamDef*);
void params_mod_split_(cseis_system::csParamDef*);
void params_mod_splitting_(cseis_system::csParamDef*);
void params_mod_stack_(cseis_system::csParamDef*);
void params_mod_statics_(cseis_system::csParamDef*);
void params_mod_sumodule_(cseis_system::csParamDef*);
void params_mod_test_(cseis_system::csParamDef*);
void params_mod_test_multi_ensemble_(cseis_system::csParamDef*);
void params_mod_test_multi_fixed_(cseis_system::csParamDef*);
void params_mod_time_slice_(cseis_system::csParamDef*);
void params_mod_time_stretch_(cseis_system::csParamDef*);
void params_mod_trc_add_ens_(cseis_system::csParamDef*);
void params_mod_trc_interpol_(cseis_system::csParamDef*);
void params_mod_trc_math_(cseis_system::csParamDef*);
void params_mod_trc_math_ens_(cseis_system::csParamDef*);
void params_mod_trc_print_(cseis_system::csParamDef*);
void params_mod_trc_split_(cseis_system::csParamDef*);

void (*METHODS_PARAM[N_METHODS])(cseis_system::csParamDef*) = {
  params_mod_attribute_,
  params_mod_bin_,
  params_mod_cmp_,
  params_mod_convolution_,
  params_mod_designature_,
  params_mod_despike_,
  params_mod_else_,
  params_mod_elseif_,
  params_mod_endif_,
  params_mod_endsplit_,
  params_mod_ens_define_,
  params_mod_fft_,
  params_mod_gain_,
  params_mod_geotools_,
  params_mod_hdr_del_,
  params_mod_hdr_math_,
  params_mod_hdr_print_,
  params_mod_hdr_set_,
  params_mod_if_,
  params_mod_image_,
  params_mod_input_,
  params_mod_input_ascii_,
  params_mod_input_create_,
  params_mod_input_rsf_,
  params_mod_input_segd_,
  params_mod_input_segy_,
  params_mod_input_sinewave_,
  params_mod_kill_,
  params_mod_lmo_,
  params_mod_mirror_,
  params_mod_mute_,
  params_mod_nmo_,
  params_mod_output_,
  params_mod_output_rsf_,
  params_mod_output_segy_,
  params_mod_p190_,
  params_mod_read_ascii_,
  params_mod_resample_,
  params_mod_rms_,
  params_mod_scaling_,
  params_mod_select_,
  params_mod_select_time_,
  params_mod_statics_,
  params_mod_test_,
  params_mod_time_stretch_,
  params_mod_trc_math_,
  params_mod_trc_print_,
  params_mod_beam_forming_,
  params_mod_ccp_,
  params_mod_concatenate_,
  params_mod_correlation_,
  params_mod_debias_,
  params_mod_fft_2d_,
  params_mod_filter_,
  params_mod_fxdecon_,
  params_mod_hdr_math_ens_,
  params_mod_histogram_,
  params_mod_hodogram_,
  params_mod_kill_ens_,
  params_mod_off2angle_,
  params_mod_orient_,
  params_mod_orient_convert_,
  params_mod_overlap_,
  params_mod_picking_,
  params_mod_poscalc_,
  params_mod_pz_sum_,
  params_mod_ray2d_,
  params_mod_repeat_,
  params_mod_resequence_,
  params_mod_rotate_,
  params_mod_semblance_,
  params_mod_sort_,
  params_mod_split_,
  params_mod_splitting_,
  params_mod_stack_,
  params_mod_sumodule_,
  params_mod_test_multi_ensemble_,
  params_mod_test_multi_fixed_,
  params_mod_time_slice_,
  params_mod_trc_add_ens_,
  params_mod_trc_interpol_,
  params_mod_trc_math_ens_,
  params_mod_trc_split_
};

#endif
