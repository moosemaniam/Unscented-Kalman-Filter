#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include "tools.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

#define DEBUG_PRINTS
/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.57;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  n_x_ = x_.size();

  n_aug_ = n_x_ + 2;

  lambda_ = 3 - n_aug_; /* From the lecture, this was the best lambda */

  n_sigma_ = 2 * n_aug_ + 1;

  Xsig_pred_ = MatrixXd(n_x_,n_sigma_);

  weights_ = VectorXd(n_sigma_);

  M_radar = MatrixXd(3,3);
  M_radar << std_radr_*std_radr_,0,0,
          0,std_radphi_*std_radphi_,0,
          0,0,std_radrd_*std_radrd_;

  M_lidar = MatrixXd(2,2);

  M_lidar << std_laspx_ * std_laspx_,0,
              0,std_laspy_*std_laspy_;

  is_initialized_=false;
#ifdef DEBUG_PRINTS
  printf("Init done\n");
#endif

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage measurement_pack) {

  /* ************************************************************ */
  /*   Initialization */
  /* ************************************************************ */
  if(!is_initialized_)
  {
    /* Set P as a identity matrix */
    P_ << 1,0,0,0,0,
       0,1,0,0,0,
       0,0,1,0,0,
       0,0,0,1,0,
       0,0,0,0,1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
        Convert radar from polar to cartesian coordinates and initialize state.
        */
      double rho = measurement_pack.raw_measurements_[0];
      double phi = measurement_pack.raw_measurements_[1];
      double rho_d= measurement_pack.raw_measurements_[2];

      /* x is projection of the rho vector on x axis */
      double x =  rho * cos(phi);
      /* x is projection of the rho vector on y axis */
      double y =  rho * sin(phi);

      /* Similiarly, velocity projections on x and y axis */
      double vx = rho_d * cos(phi);
      double vy = rho_d * sin(phi);
      double v = sqrt(vx*vx + vy*vy);

      x_ << x,y,v,0,0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
        Initialize state.
        */
      x_ << measurement_pack.raw_measurements_[0],measurement_pack.raw_measurements_[1],0,0,0;
    }

    if(fabs(x_[0]) < MIN_VAL)
      x_[0] = MIN_VAL;

    if(fabs(x_[1]) < MIN_VAL)
      x_[1] = MIN_VAL;

    weights_(0) = lambda_/(n_aug_ + lambda_);

    for(int i=1;i < weights_.size();i++)
    {
      weights_(i) = 0.5/(n_aug_ + lambda_);
    }
    /* First time stamp measurement */
    time_us_ = measurement_pack.timestamp_;
    is_initialized_ = true;

    /* After init, come out. Update and prediction will follow */
    /* from the next input onwards */
    return;
  }

  double dt =  measurement_pack.timestamp_ - time_us_;
  time_us_ = measurement_pack.timestamp_;

  /* Convert to micro-seconds */
  dt = dt/1000000.0;

  Prediction(dt);



  if(measurement_pack.sensor_type_ == MeasurementPackage::RADAR &&
      use_radar_) {
    UpdateRadar(measurement_pack);

  }
  if(measurement_pack.sensor_type_ == MeasurementPackage::LASER &&
      use_laser_) {
    UpdateLidar(measurement_pack);

  }
  //exit(0);
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  Tools math_tools;
  double dt2 = delta_t * delta_t;
  double sqrt_lambda_plus_n_aug = sqrt(lambda_ + n_aug_);
  VectorXd sqrt_lambda_plus_n_aug_vec;


  VectorXd x_aug = VectorXd(n_aug_); /* Augmented mean vector */
  MatrixXd Xsig_aug = MatrixXd(n_aug_,n_sigma_); /* Sigma point matrix */
  MatrixXd P_aug = MatrixXd(n_aug_,n_aug_); /* Augmented state covariance matrix */
  /* Fill x aug */
  Xsig_aug.fill(0.0);
  x_aug.fill(0.0);
  x_aug.head(n_x_) = x_;

  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_,n_x_) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_* std_yawdd_;

  /* Find Square root of the P matrix */
  MatrixXd L = P_aug.llt().matrixL();

  Xsig_aug.col(0) = x_aug;

  for(int i=0;i<n_aug_;i++){
    sqrt_lambda_plus_n_aug_vec = sqrt_lambda_plus_n_aug * L.col(i);
    Xsig_aug.col(i+1) = x_aug + sqrt_lambda_plus_n_aug_vec;
    Xsig_aug.col(i+n_aug_+1) = x_aug - sqrt_lambda_plus_n_aug_vec;

  }

  /* Sigma points predictions */
  for(int i =0; i< n_sigma_;i++)
  {
    double p_x= Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v= Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yaw_d = Xsig_aug(4,i);

    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);
    double update_angle = yaw + yaw_d*delta_t;

    double  px_predicted,py_predicted;

    double vx_yawd,px_update,py_update;


    if(fabs(yaw_d) > MIN_VAL)
    { double v_by_yawd = v / yaw_d;
      px_predicted = p_x + v_by_yawd * (sin(update_angle) - sin(yaw));
      py_predicted = p_y + v_by_yawd * (cos(yaw) - cos(update_angle));

    }
    else {
      /* Incase of low/zero yaw_d value */
      px_predicted= p_x + v  * delta_t * cos(yaw);
      py_predicted= p_y + v  * delta_t * sin(yaw);
    }

    double v_p = v;
    double yaw_p,yawd_p;
    /*Position Prediction */

    yawd_p = yaw_d;
    yaw_p = update_angle;
    /* add noise to position prediction */
    px_predicted += 0.5*nu_a*dt2*cos(yaw);
    py_predicted += 0.5*nu_a*dt2*sin(yaw);


    /* Velocity and yaw predictions */
    v_p +=  nu_a * delta_t;
    yaw_p += 0.5*nu_yawdd*dt2;
    yawd_p += nu_yawdd*delta_t;

    /* Update predicted sigma points */
    Xsig_pred_(0,i) = px_predicted;
    Xsig_pred_(1,i) = py_predicted;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

  x_ =  Xsig_pred_ * weights_;

  P_.fill(0.0);

  for(int i=0;i < n_sigma_;i++){
    VectorXd x_d = Xsig_pred_.col(i) - x_;
    x_d(3) =  math_tools.normalize_to_pi(x_d(3));
    P_ = P_ + weights_(i) * x_d * x_d.transpose();
  }
}
/* Common update function for both lidar and radar */
void UKF::UKFCommonUpdate(MeasurementPackage meas_package,
    MatrixXd Sig,
    int measurement_size)
{
  VectorXd pred = VectorXd(measurement_size);
  pred = Sig * weights_;
  Tools math_tools;

  /* measurement covariance matrix*/
  MatrixXd S = MatrixXd(measurement_size,measurement_size);
  S.fill(0.0);
  for(int i =0;i < n_sigma_;i++)
  {
    VectorXd diff = Sig.col(i) - pred;
    /* normalize the angle to -pi to pi*/
    diff(1) = math_tools.normalize_to_pi(diff(1));
    S = S + weights_(i) * diff * diff.transpose();
  }

  MatrixXd M_noise = MatrixXd(measurement_size,measurement_size);

  if(meas_package.sensor_type_ ==  MeasurementPackage::RADAR)
  {
    M_noise = M_radar;
  }
  else
  {
    M_noise = M_lidar;
  }

  /* Add measurement noise */
  S = S + M_noise;
  /* Cross correlation matrix */
  MatrixXd CC_M = MatrixXd(n_x_,measurement_size);
  CC_M.fill(0.0);

  for(int i =0;i < n_sigma_;i++)
  {
    VectorXd diff = Sig.col(i) - pred;
    if(meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
      diff(1) = math_tools.normalize_to_pi(diff(1));
    }
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    x_diff(3) = math_tools.normalize_to_pi(x_diff(3));
    CC_M = CC_M + weights_(i) * x_diff * diff.transpose();
  }

  VectorXd z = meas_package.raw_measurements_;
  MatrixXd KalmanGain =  CC_M * S.inverse();

  VectorXd diff = z - pred;

  /* Normalize angle for radar */
if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
{
    diff(1) = math_tools.normalize_to_pi(diff(1));
  }

  x_ = x_ + KalmanGain * diff;
  P_ = P_ - KalmanGain * S * KalmanGain.transpose();

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {
    nis_radar = z.transpose() * S.inverse() * z;
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER)
  {
	  nis_lidar= z.transpose() * S.inverse() * z;
  }
}
/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  const int measurement_size=2;
  MatrixXd XsigLidar = Xsig_pred_.block(0,0,measurement_size,n_sigma_);
  UKFCommonUpdate(meas_package,XsigLidar,measurement_size);
}



/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  const int measurement_size=3;
#if 1

  MatrixXd XsigRadar= MatrixXd(measurement_size,n_sigma_);

  XsigRadar.fill(0.0);
  for(int i=0;i < n_sigma_;i++){
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    /* horizontal component of velocity */
    double yaw = Xsig_pred_(3,i);
    double v_x = v*cos(yaw);
    /* vertical component of velocity */
    double v_y = v*sin(yaw);

    /* r */
    XsigRadar(0,i) =sqrt(p_x*p_x + p_y*p_y);
    /* phi */
    XsigRadar(1,i) =atan2(p_y,p_x);
    /* dot */
    XsigRadar(2,i) =(p_x*v_x + p_y*v_y)/XsigRadar(0,i);
  }
  UKFCommonUpdate(meas_package,XsigRadar,measurement_size);
#endif
}
