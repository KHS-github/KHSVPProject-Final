using ROOT::Math;

const double boxL = 1.0;

class Particle
{
public:
 Particle();
 ~Particle();
private:
 double _mass;
 ROOT::Math::XYZVector _position;
 ROOT::Math::XYZVector _momentum;
public:
 void move(double delta_time);
 bool collision_detector(XYZVector& delta_pos);
};

Particle::Particle(double mass, XYZVector position, XYZVector momentum) : _mass(mass), _position(position), _momentum(momentum)
{
}

Particle::~Particle()
{
}

void Particle::move(double delta_time)
{
 auto delta_pos = _momentum / _mass * delta_time;
 auto fut = _position + delta_pos;

 double collision_plane[] = {0.0, 0.0, 0.0, 0.0};

 if(delta_pos.x() < 0.0){
  collision_plane[0] = 1.0;
  collision_plane[3] = 0.0;
 }
 else if(delta_pos.x() > boxL){
  collision_plane[0] = -1.0;
  collision_plane[3] = boxL;
 }
 else if(delta_pos.y() < 0.0){
  collision_plane[1] = 1.0;
  collision_plane[3] = 0.0;
 }
 else if(delta_pos.y() > boxL){
  collision_plane[1] = -1.0;
  collision_plane[3] = boxL;
 }
 else if(delta_pos.z() < 0.0){
  collision_plane[2] = 1.0;
  collision_plane[3] = 0.0;
 }
 else if(delta_pos.z() > boxL){
  collision_plane[2] = -1.0;
  collision_plane[3] = 0.0;
 }
 else{
  _position += delta_pos;
  return;
 }

 //get collision time
 double _4d_pos[] = {_position.x(), _position.y(), _position.z(), 1.0};
 double _4d_vel[] = {_momentum.x()/_mass, _momentum.y()/_mass, _momentum.z()/_mass, 0.0};
 double delta_t_star = -getDot4D(collision_plane, _4d_pos) / getDot4D(collision_plane, _4d_vel);
 
 //get collided plane point
 XYZVector new_delta = delta_t_star * _momentum / _mass;
 

 //get reversed momentum
 XYZVector plane_normal(collision_plane[0], collision_plane[1], collision_plane[2]);
 double normal_contrib = plane_normal.Dot(_momentum);
 _momentum += -normal_contrib * plane_normal;

 //update position
 double remain_time = delta_t - delta_t_star; 
 new_delta += _momentum / _mass * remain_time;
 _position += new_delta;
}

double getDot4D(double* arr1, double* arr2)
{
 double ret = 0.0;
 for(int i=0;i<4;i++){
  ret += arr1[i] * arr2[i];
 }
 return ret;
}

void ergodicity()
{
 const double ESTIMATED_TOTAL_TIME = 1000;
 const double DELTA_TIME = 0.00001;
 const int TOTAL_STEPS = (int)(ESTIMATED_TOTAL_TIME / DELTA_TIME) + 1;
 list<Particle*> particles;
 
 //create particles and set initial onditions to them
 const int PARTICLE_NUM = 1;
 for(int i=0;i<PARTICLE_NUM;i++){
  particles.push_back(1.0, XYZVector(), XYZVector());
 }


 
 for(int i=0;i<TOTAL_STEPS;i++){
  for(auto particle : particles)
  {
   particle.move(DETLA_TIME);
  }
 }
}
