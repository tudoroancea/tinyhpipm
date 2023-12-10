#ifdef __cplusplus
extern "C" {
#endif

void m_fact_solve_kkt_step_hard_ocp_qp(struct d_ocp_qp* d_qp, struct s_ocp_qp* s_qp, struct m_ipm_hard_ocp_qp_workspace* ws);
void m_solve_kkt_step_hard_ocp_qp(struct d_ocp_qp* d_qp, struct s_ocp_qp* s_qp, struct m_ipm_hard_ocp_qp_workspace* ws);

#ifdef __cplusplus
} /* extern "C" */
#endif
