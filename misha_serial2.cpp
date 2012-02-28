#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"

//================================================================================
//=========================== Global Constants ===================================
//================================================================================

//Particle Radius
double radius = cutoff/2;

//Token types
const int R = 1;
const int L = -1;

//================================================================================
//========================== Structures ==========================================
//================================================================================

//Collision Token
typedef struct {
  double position;
  int type;
  int particle_id;
} token;

//Active collisions
//List of [particle_id1 particle_id2] pairs
typedef struct {
  int id1;
  int id2;
} collision;

//Holds the information for a single partition
typedef struct {
  //Particles
  particle_t* particles;
  int num_particles;
  int max_particles;

  //Active IDs
  int* is_id_active;
  int* free_ids;

  //Collision Tokens
  token* xtokens;
  token* ytokens;

  //Active Collisions
  collision* active_collisions;
  int num_active_collisions;

  //Collision table
  int* collision_table;
} partition;

//================================================================================
//========================== Collision Detector Fields  ==========================
//================================================================================

//Particles
particle_t* particles;
int num_particles;
int max_particles;

//Collision Tokens
token* xtokens;
token* ytokens;

//Active Collisions
collision* active_collisions;
int num_active_collisions;

//Collision Table
int* collision_table;

//Active IDs
int* is_id_active;
int* free_ids;
int num_used_ids;

//================================================================================
//========================== Triangular Matrix Utilities =========================
//================================================================================

//Compute the number of elements in an upper triangular matrix excluding the diagonal
inline int num_triangle_elements(int n){
  return (n-1)*n/2;
}

//Compute the idx of the element at (row,col) in an upper triangular matrix excluding diagonal
inline int triangle_idx(int n, int row, int col){
  return n*row - (row+1)*(row+2)/2 + col;
}

//================================================================================
//========================= Partition Manager ====================================
//================================================================================

//Set partition p as the active partition
void set_active_partition(partition* p){
  //Particles
  particles = p->particles;
  num_particles = p->num_particles;
  max_particles = p->max_particles;

  //Collision Tokens
  xtokens = p->xtokens;
  ytokens = p->ytokens;

  //Active Collisions
  active_collisions = p->active_collisions;
  num_active_collisions = p->num_active_collisions;

  //Collision Table
  collision_table = p->collision_table;

  //Active IDs
  is_id_active = p->is_id_active;
  free_ids = p->free_ids;
}

//Create a new partition
partition* alloc_partition(int max_particles){
  //Allocate partition structure
  partition* p = (partition*)malloc(sizeof(partition));

  //Allocate space for particles
  p->max_particles = max_particles;
  p->num_particles = 0;
  p->particles = (particle_t*)malloc(max_particles * sizeof(particle_t));

  //Allocate buffers
  p->xtokens = (token*)malloc(2 * max_particles * sizeof(token));
  p->ytokens = (token*)malloc(2 * max_particles * sizeof(token));

  //Allocate active collision list  
  p->num_active_collisions = 0;
  p->active_collisions = (collision*)malloc(max_particles * 2 * sizeof(collision));

  //Allocate collision table
  int size = num_triangle_elements(max_particles);
  p->collision_table = (int*)malloc(size * sizeof(int));
  for(int i=0; i<size; i++)
    p->collision_table[i] = 0;

  //Active and Free IDs
  p->is_id_active = (int*)malloc(max_particles * sizeof(int));
  p->free_ids = (int*)malloc(max_particles * sizeof(int));
  for(int i=0; i<max_particles; i++){
    p->is_id_active[i] = 0;
    p->free_ids[i] = i;
  }

  //Return new partition
  return p;
}

//================================================================================
//======================= Collision Detector =====================================
//================================================================================

//Register an active collision
void register_active_collision(int id1, int id2){
  collision* c = &active_collisions[num_active_collisions++];
  c->id1 = id1;
  c->id2 = id2;
}

//Mark id1 and id2 as intersecting
void mark_intersection(int id1, int id2){
  int min_id = min(id1, id2);
  int max_id = max(id1, id2);
  
  int idx = triangle_idx(max_particles, min_id, max_id);
  int num_intersections = collision_table[idx];
  collision_table[idx] = num_intersections + 1;

  if(num_intersections == 1)
    register_active_collision(min_id, max_id);
}

//Unmark id1 and id2 as intersecting
void unmark_intersection(int id1, int id2){
  int min_id = min(id1, id2);
  int max_id = max(id1, id2);
  
  int idx = triangle_idx(max_particles, min_id, max_id);
  collision_table[idx]--;
}

//Notify collision detector that token t1 and t2 has been swapped
void swap(token t1, token t2){
  if(t1.type == R && t2.type == L)
    mark_intersection(t1.particle_id, t2.particle_id);
  else if(t1.type == L && t2.type == R)
    unmark_intersection(t1.particle_id, t2.particle_id);
}

//Sinking sort with collision detection logic
//Sinks element i down to its rightful position.
//Assumes list up to i is sorted.
inline void sweep_down(token* tokens, int i){
  while(i>0) {
    token t1 = tokens[i-1];
    token t2 = tokens[i];
    
    if(!is_id_active[t1.particle_id]){
      tokens[i] = t1;
      tokens[i-1] = t2;
      i--;
    } else if(t2.position < t1.position){
      //Swap
      tokens[i] = t1;
      tokens[i-1] = t2;
      //Notify swap
      swap(t1, t2);
      i--;
    } else{
      //Found rightful position
      break;
    }
  }
}

//Sweep through tokens from i = 1 ... num_tokens
//Sinking sort elements
void sweep(token* tokens){
  for(int i=1; i<2 * num_particles; i++)
    if(is_id_active[tokens[i].particle_id])
       sweep_down(tokens, i);
}

//Prune the active collisions
//Keep only the ones where num intersections is >= 2
void prune(){
  collision* dest = active_collisions;

  int collision_length = num_active_collisions;  
  for(int i=0; i<collision_length; i++){
    collision* c = &active_collisions[i];

    //Check whether ids are active
    int is_active = is_id_active[c->id1] && is_id_active[c->id2];
    
    //Check whether collision is active
    int idx = triangle_idx(max_particles, c->id1, c->id2);
    int num_intersections = collision_table[idx];

    //If the collision is active then copy to dest
    if(is_active && num_intersections == 2){
      if(dest != c)
        *dest = *c;
      dest++;
    }else{
      //Otherwise, discount this collision
      num_active_collisions--;
    }
  }
}

//Update the token information in the collision detector
void update_tokens(){
  //Update x tokens
  for(int i=0; i<2 * num_particles; i++){
    token* t = &xtokens[i];
    if(is_id_active[t->particle_id]){
      particle_t p = particles[t->particle_id];
      t->position = p.x + (t->type * radius);
    }
  }

  //Update y tokens
  for(int i=0; i<2 * num_particles; i++){
    token* t = &ytokens[i];
    if(is_id_active[t->particle_id]){
      particle_t p = particles[t->particle_id];
      t->position = p.y + (t->type * radius);
    }
  }
}

//Full Collision Detection Sweep
void sweep_and_prune(){
  update_tokens();
  sweep(xtokens);
  sweep(ytokens);
  prune();
}

//================================================================================
//====================== Collision Detector Interface ============================
//================================================================================

int add_particle(particle_t* p){
  //Get id
  int id = free_ids[num_used_ids++];
  //Set active
  is_id_active[id] = 1;
  //Copy data
  particles[id] = *p;
  
  //Create x tokens
  token* L_xtoken = &xtokens[2 * num_particles];
  L_xtoken->type = L;
  L_xtoken->particle_id = id;
  
  token* R_xtoken = &xtokens[2 * num_particles + 1];
  R_xtoken->type = R;
  R_xtoken->particle_id = id;

  //Create y tokens
  token* L_ytoken = &ytokens[2 * num_particles];
  L_ytoken->type = L;
  L_ytoken->particle_id = id;
  
  token* R_ytoken = &ytokens[2 * num_particles + 1];
  R_ytoken->type = R;
  R_ytoken->particle_id = id;

  //Increment num_particles
  num_particles++;

  //Clear Collision Table
  for(int i=0; i<num_particles; i++){
    int idx = triangle_idx(max_particles, min(i,id), max(i,id));
    collision_table[idx] = 0;
  }
}

void remove_particle(int id){
  //Set inactive
  is_id_active[id] = 0;
  //Return id
  num_used_ids--;
  free_ids[num_used_ids] = id;
}

void set_state(int id, double x, double y, double vx, double vy){
  if(!is_id_active[id]){
    printf("Particle %d is not active.\n", id);
    exit(-1);
  }
  particles[id].x = x;
  particles[id].y = y;
  particles[id].vx = vx;
  particles[id].vy = vy;
}

particle_t* get_particle(int id){
  if(!is_id_active[id]){
    printf("Particle %d is not active.\n", id);
    exit(-1);
  }
  return &particles[id];
}

//================================================================================
//==================== Sorted Initialization =====================================
//================================================================================

//  //Simple insertion sort
//  //Sort ps according to p.x or p.y depending on yaxis
//  //Returns sorted ordering in order
//  int sort(particle_t* ps, int* order, int n, int yaxis){
//    for(int i=1; i<n; i++){
//      int idx = order[i];
//      double xi = yaxis? ps[idx].y : ps[idx].x;
//      
//      int j = i-1;
//      double xj = yaxis? ps[order[j]].y : ps[order[j]].x;    
//      while(j>=0 && xj > xi){
//        order[j+1] = order[j];
//        j--;
//        if(j>=0)
//           xj = yaxis? ps[order[j]].y : ps[order[j]].x;
//      }
//      order[j+1] = idx;
//    }
//  }
//  
//  //Setup the collision detector
//  void setup_collision_detector(particle_t* ps, int n){
//    for(int i=0; i<n; i++)
//      add_particle(&ps[i]);
//  
//    /*
//    
//    //Register particles
//    num_particles = n;
//    for(int i=0; i<n; i++)
//      particles[i] = ps[i];
//  
//    //Sorted Initialization
//    int* order = (int*)malloc(num_particles * sizeof(int));
//    for(int i=0; i<num_particles; i++)
//      order[i] = i;
//  
//    //Create x tokens
//    sort(particles, order, num_particles, 0);
//    for(int i=0, j=0; j<num_particles; i+=2, j++){
//      //L xtoken
//      token* L_xtoken = &xtokens[i];
//      L_xtoken->type = L;
//      L_xtoken->particle_id = order[j];
//      
//      //R xtoken
//      token* R_xtoken = &xtokens[i+1];
//      R_xtoken->type = R;
//      R_xtoken->particle_id = order[j];
//    }
//  
//    //Create y tokens
//    sort(particles, order, num_particles, 1);
//    for(int i=0, j=0; j<num_particles; i+=2, j++){
//      //L ytoken
//      token* L_ytoken = &ytokens[i];
//      L_ytoken->type = L;
//      L_ytoken->particle_id = order[j];
//      
//      //R ytoken
//      token* R_ytoken = &ytokens[i+1];
//      R_ytoken->type = R;
//      R_ytoken->particle_id = order[j];
//    }
//    free(order);
//  
//    */
//  }

//================================================================================
//================== Physics Calculations ========================================
//================================================================================

//Apply a force to p1, and an equal and opposing force (ala Newton's 3rd law) to p2
void apply_pairwise_force(particle_t* p1, particle_t* p2) {
  double dx = p2->x - p1->x;
  double dy = p2->y - p1->y;
  double r2 = dx * dx + dy * dy;
  if(r2 < cutoff*cutoff) {
      
    //Limit the maximum force
    r2 = fmax( r2, min_r*min_r );
    double r = sqrt( r2 );

    //Repulsive Force
    double coef = ( 1 - cutoff / r ) / r2 / mass;

    //Propel both particles
    p1->ax += coef * dx;
    p1->ay += coef * dy;
    p2->ax -= coef * dx;
    p2->ay -= coef * dy;
  }
}

void update_particles(){
  //Calculate active collisions
  sweep_and_prune();

  //Reset acceleration
  for(int i=0; i<num_particles; i++)
    if(is_id_active[i])
      particles[i].ax = particles[i].ay = 0;

  //Accumulate acceleration
  for(int i=0; i<num_active_collisions; i++){
    collision c = active_collisions[i];
    apply_pairwise_force(&particles[c.id1], &particles[c.id2]);
  }

  //Move Particles
  for(int i=0; i<num_particles; i++)
    if(is_id_active[i])
      move(particles[i]);
}

void run_simulation(particle_t* ps, int n, FILE* fsave){
  //Create partition and set active
  partition* p = alloc_partition(n);
  set_active_partition(p);

  //Add all particles
  for(int i=0; i<n; i++)
    add_particle(&ps[i]);

  //For each step
  for(int step = 0; step < NSTEPS; step++ ){
    update_particles();
    
    //=== Save state to file ===
    if( fsave && (step%SAVEFREQ) == 0 )
      save( fsave, n, particles );
  }
}

//================================================================================
//==================== Main Driver ===============================================
//================================================================================

int main( int argc, char **argv )
{    
    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        return 0;
    }
    
    int n = read_int( argc, argv, "-n", 1000 );

    char *savename = read_string( argc, argv, "-o", NULL );
    
    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    set_size( n );
    init_particles( n, particles );
    
    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
    run_simulation(particles, n, fsave);
    simulation_time = read_timer( ) - simulation_time;
    
    printf( "n = %d, simulation time = %g seconds\n", n, simulation_time );
    
    free( particles );
    if( fsave )
        fclose( fsave );
    
    return 0;
}
