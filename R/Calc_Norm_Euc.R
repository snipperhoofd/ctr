#' Calculate Euclidean Distance
#'
#' \deqn{$d(q,p)=\sqrt{(q_{1}-p_{1})^2+...+(q_{n}-p_{n})^2}$}
#'
#' {\displaystyle {\begin{aligned}\mathrm {d} (\mathbf {p} ,\mathbf {q} )=\mathrm {d} (\mathbf {q} ,\mathbf {p} )&={\sqrt {(q_{1}-p_{1})^{2}+(q_{2}-p_{2})^{2}+\cdots +(q_{n}-p_{n})^{2}}}\\[8pt]&={\sqrt {\sum _{i=1}^{n}(q_{i}-p_{i})^{2}}}.\end{aligned}}} \begin{align}\mathrm{d}(\mathbf{p},\mathbf{q}) = \mathrm{d}(\mathbf{q},\mathbf{p}) & = \sqrt{(q_1-p_1)^2 + (q_2-p_2)^2 + \cdots + (q_n-p_n)^2} \\[8pt]
#' & = \sqrt{\sum_{i=1}^n (q_i-p_i)^2}.\end{align}
#'
#' @export
#' @param vector1,vector2  Two vectors for which the euclidean disance is to be calculated
#'
#' @return The Euclidean Distance between two two vectors.
#' @examples Euclidean_Distance<-Calc_Norm_Euc(as.numeric(RNAseq_Annotation_Matrix_no_sd_of_zero[position_of_A,RS:RE]),as.numeric(RNAseq_Annotation_Matrix_no_sd_of_zero[position_of_B,RS:RE]))

Calc_Norm_Euc <- function(vector1, vector2){
  Norm_Euc_Distance<-sqrt(sum((vector1-vector2)*(vector1-vector2)))
  return(Norm_Euc_Distance)
}
