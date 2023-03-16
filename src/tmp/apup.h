arma::vec a_up1(const arma::vec & p, const int & ind1, const int & ind2);

arma::vec a_up3(const arma::vec & p, const int & Lind1, const int & Lind2, const int & Rind1, const int & Rind2);

double sumxlogy(const arma::vec & x, const arma::vec & y);

double p_up(arma::vec & p, const arma::vec & d, const double & alpha0);

arma::umat acount(const arma::vec & L, const arma::vec & R, const arma::vec & breaks);
