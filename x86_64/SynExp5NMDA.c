/* Created by Language version: 7.7.0 */
/* NOT VECTORIZED */
#define NRN_VECTORIZED 0
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib_ansi.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__Exp5NMDA
#define _nrn_initial _nrn_initial__Exp5NMDA
#define nrn_cur _nrn_cur__Exp5NMDA
#define _nrn_current _nrn_current__Exp5NMDA
#define nrn_jacob _nrn_jacob__Exp5NMDA
#define nrn_state _nrn_state__Exp5NMDA
#define _net_receive _net_receive__Exp5NMDA 
#define rates rates__Exp5NMDA 
#define state state__Exp5NMDA 
 
#define _threadargscomma_ /**/
#define _threadargsprotocomma_ /**/
#define _threadargs_ /**/
#define _threadargsproto_ /**/
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 static double *_p; static Datum *_ppvar;
 
#define t nrn_threads->_t
#define dt nrn_threads->_dt
#define tau1 _p[0]
#define tau2_0 _p[1]
#define a2 _p[2]
#define b2 _p[3]
#define wtau2 _p[4]
#define tau3_0 _p[5]
#define a3 _p[6]
#define b3 _p[7]
#define tp _p[8]
#define d1 _p[9]
#define tau_D1 _p[10]
#define tauV _p[11]
#define gVDst _p[12]
#define gVDv0 _p[13]
#define gVI _p[14]
#define Mg _p[15]
#define K0 _p[16]
#define delta _p[17]
#define e _p[18]
#define i _p[19]
#define wf _p[20]
#define A _p[21]
#define B _p[22]
#define C _p[23]
#define gVD _p[24]
#define g _p[25]
#define factor _p[26]
#define q10_tau2 _p[27]
#define q10_tau3 _p[28]
#define tau _p[29]
#define wtau3 _p[30]
#define DA _p[31]
#define DB _p[32]
#define DC _p[33]
#define DgVD _p[34]
#define _g _p[35]
#define _tsav _p[36]
#define _nd_area  *_ppvar[0]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static double _hoc_Mgblock();
 static double _hoc_rates();
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern Prop* nrn_point_prop_;
 static int _pointtype;
 static void* _hoc_create_pnt(_ho) Object* _ho; { void* create_point_process();
 return create_point_process(_pointtype, _ho);
}
 static void _hoc_destroy_pnt();
 static double _hoc_loc_pnt(_vptr) void* _vptr; {double loc_point_process();
 return loc_point_process(_pointtype, _vptr);
}
 static double _hoc_has_loc(_vptr) void* _vptr; {double has_loc_point();
 return has_loc_point(_vptr);
}
 static double _hoc_get_loc_pnt(_vptr)void* _vptr; {
 double get_loc_point_process(); return (get_loc_point_process(_vptr));
}
 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _p = _prop->param; _ppvar = _prop->dparam;
 }
 static void _hoc_setdata(void* _vptr) { Prop* _prop;
 _prop = ((Point_process*)_vptr)->_prop;
   _setdata(_prop);
 }
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 0,0
};
 static Member_func _member_func[] = {
 "loc", _hoc_loc_pnt,
 "has_loc", _hoc_has_loc,
 "get_loc", _hoc_get_loc_pnt,
 "Mgblock", _hoc_Mgblock,
 "rates", _hoc_rates,
 0, 0
};
#define Mgblock Mgblock_Exp5NMDA
 extern double Mgblock( double );
 /* declare global and static user variables */
#define Q10 Q10_Exp5NMDA
 double Q10 = 1.52;
#define Q10_tau3 Q10_tau3_Exp5NMDA
 double Q10_tau3 = 2.65;
#define Q10_tau2 Q10_tau2_Exp5NMDA
 double Q10_tau2 = 3.68;
#define Q10_tau1 Q10_tau1_Exp5NMDA
 double Q10_tau1 = 2.2;
#define T0 T0_Exp5NMDA
 double T0 = 26;
#define T0_tau T0_tau_Exp5NMDA
 double T0_tau = 35;
#define inf inf_Exp5NMDA
 double inf = 0;
#define tau3 tau3_Exp5NMDA
 double tau3 = 0;
#define tau2 tau2_Exp5NMDA
 double tau2 = 0;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "d1", 0, 1,
 "tau_D1", 1e-09, 1e+09,
 "tauV", 1e-09, 1e+09,
 "tau1", 1e-09, 1e+09,
 "wtau2", 1e-09, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "T0_tau_Exp5NMDA", "degC",
 "T0_Exp5NMDA", "degC",
 "inf_Exp5NMDA", "uS",
 "tau2_Exp5NMDA", "ms",
 "tau3_Exp5NMDA", "ms",
 "tau1", "ms",
 "tau2_0", "ms",
 "a2", "ms",
 "b2", "1/mV",
 "wtau2", "1e-9",
 "tau3_0", "ms",
 "a3", "ms",
 "b3", "1/mV",
 "tp", "ms",
 "d1", "1",
 "tau_D1", "ms",
 "tauV", "ms",
 "gVDst", "1/mV",
 "gVDv0", "mV",
 "gVI", "uS",
 "Mg", "mM",
 "K0", "mM",
 "delta", "1",
 "e", "mV",
 "gVD", "uS",
 "i", "nA",
 0,0
};
 static double A0 = 0;
 static double B0 = 0;
 static double C0 = 0;
 static double delta_t = 0.01;
 static double gVD0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "Q10_tau1_Exp5NMDA", &Q10_tau1_Exp5NMDA,
 "Q10_tau2_Exp5NMDA", &Q10_tau2_Exp5NMDA,
 "Q10_tau3_Exp5NMDA", &Q10_tau3_Exp5NMDA,
 "T0_tau_Exp5NMDA", &T0_tau_Exp5NMDA,
 "Q10_Exp5NMDA", &Q10_Exp5NMDA,
 "T0_Exp5NMDA", &T0_Exp5NMDA,
 "inf_Exp5NMDA", &inf_Exp5NMDA,
 "tau2_Exp5NMDA", &tau2_Exp5NMDA,
 "tau3_Exp5NMDA", &tau3_Exp5NMDA,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void nrn_cur(_NrnThread*, _Memb_list*, int);
static void  nrn_jacob(_NrnThread*, _Memb_list*, int);
 static void _hoc_destroy_pnt(_vptr) void* _vptr; {
   destroy_point_process(_vptr);
}
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(_NrnThread*, _Memb_list*, int);
static void _ode_matsol(_NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[2]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"Exp5NMDA",
 "tau1",
 "tau2_0",
 "a2",
 "b2",
 "wtau2",
 "tau3_0",
 "a3",
 "b3",
 "tp",
 "d1",
 "tau_D1",
 "tauV",
 "gVDst",
 "gVDv0",
 "gVI",
 "Mg",
 "K0",
 "delta",
 "e",
 0,
 "i",
 "wf",
 0,
 "A",
 "B",
 "C",
 "gVD",
 0,
 0};
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
  if (nrn_point_prop_) {
	_prop->_alloc_seq = nrn_point_prop_->_alloc_seq;
	_p = nrn_point_prop_->param;
	_ppvar = nrn_point_prop_->dparam;
 }else{
 	_p = nrn_prop_data_alloc(_mechtype, 37, _prop);
 	/*initialize range parameters*/
 	tau1 = 1.69;
 	tau2_0 = 3.97;
 	a2 = 0.7;
 	b2 = 0.0243;
 	wtau2 = 0.65;
 	tau3_0 = 41.62;
 	a3 = 34.69;
 	b3 = 0.01;
 	tp = 30;
 	d1 = 0.2;
 	tau_D1 = 2500;
 	tauV = 7;
 	gVDst = 0.007;
 	gVDv0 = -100;
 	gVI = 1;
 	Mg = 1;
 	K0 = 4.1;
 	delta = 0.8;
 	e = -0.7;
  }
 	_prop->param = _p;
 	_prop->param_size = 37;
  if (!nrn_point_prop_) {
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 3, _prop);
  }
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _net_receive(Point_process*, double*, double);
 static void _net_init(Point_process*, double*, double);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _SynExp5NMDA_reg() {
	int _vectorized = 0;
  _initlists();
 	_pointtype = point_register_mech(_mechanism,
	 nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init,
	 hoc_nrnpointerindex, 0,
	 _hoc_create_pnt, _hoc_destroy_pnt, _member_func);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 37, 3);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
  hoc_register_dparam_semantics(_mechtype, 2, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 pnt_receive[_mechtype] = _net_receive;
 pnt_receive_init[_mechtype] = _net_init;
 pnt_receive_size[_mechtype] = 3;
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 Exp5NMDA /Users/hadasmanor/Documents/GitHubProjects/OOP_SimulationGUI/x86_64/SynExp5NMDA.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double T = 273.16;
 static double F = 9.648e4;
 static double R = 8.315;
 static double z = 2;
static int _reset;
static char *modelname = "Triple-exp model of NMDAR has (HH-type gating) (temp. sensitivity) (voltage-dependent time constants) (desensitization)";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int rates(double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static double *_temp1;
 static int _slist1[4], _dlist1[4];
 static int state(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 () {_reset=0;
 {
   rates ( _threadargscomma_ v ) ;
   DA = - A / tau1 ;
   DB = - B / tau2 ;
   DC = - C / tau3 ;
   DgVD = ( ( wtau3 * C + wtau2 * B ) / wf ) * ( inf - gVD ) / tau ;
   }
 return _reset;
}
 static int _ode_matsol1 () {
 rates ( _threadargscomma_ v ) ;
 DA = DA  / (1. - dt*( ( - 1.0 ) / tau1 )) ;
 DB = DB  / (1. - dt*( ( - 1.0 ) / tau2 )) ;
 DC = DC  / (1. - dt*( ( - 1.0 ) / tau3 )) ;
 DgVD = DgVD  / (1. - dt*( ( ( ( ( wtau3 * C + wtau2 * B ) / wf ) )*( ( ( - 1.0 ) ) ) ) / tau )) ;
  return 0;
}
 /*END CVODE*/
 
static int state () {_reset=0;
 {
   rates ( _threadargscomma_ v ) ;
   DA = - A / tau1 ;
   DB = - B / tau2 ;
   DC = - C / tau3 ;
   DgVD = ( ( wtau3 * C + wtau2 * B ) / wf ) * ( inf - gVD ) / tau ;
   }
 return _reset;}
 
static void _net_receive (_pnt, _args, _lflag) Point_process* _pnt; double* _args; double _lflag; 
{    _p = _pnt->_prop->param; _ppvar = _pnt->_prop->dparam;
  if (_tsav > t){ extern char* hoc_object_name(); hoc_execerror(hoc_object_name(_pnt->ob), ":Event arrived out of order. Must call ParallelContext.set_maxstep AFTER assigning minimum NetCon.delay");}
 _tsav = t; {
   _args[1] = 1.0 - ( 1.0 - _args[1] ) * exp ( - ( t - _args[2] ) / tau_D1 ) ;
   _args[2] = t ;
   wf = _args[0] * factor * _args[1] ;
     if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for general derivimplicit and KINETIC case */
    int __i, __neq = 4;
    double __state = A;
    double __primary_delta = (A + wf) - __state;
    double __dtsav = dt;
    for (__i = 0; __i < __neq; ++__i) {
      _p[_dlist1[__i]] = 0.0;
    }
    _p[_dlist1[0]] = __primary_delta;
    dt *= 0.5;
    v = NODEV(_pnt->node);
    _ode_matsol_instance1(_threadargs_);
    dt = __dtsav;
    for (__i = 0; __i < __neq; ++__i) {
      _p[_slist1[__i]] += _p[_dlist1[__i]];
    }
  } else {
 A = A + wf ;
     }
   if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for general derivimplicit and KINETIC case */
    int __i, __neq = 4;
    double __state = B;
    double __primary_delta = (B + wf) - __state;
    double __dtsav = dt;
    for (__i = 0; __i < __neq; ++__i) {
      _p[_dlist1[__i]] = 0.0;
    }
    _p[_dlist1[1]] = __primary_delta;
    dt *= 0.5;
    v = NODEV(_pnt->node);
    _ode_matsol_instance1(_threadargs_);
    dt = __dtsav;
    for (__i = 0; __i < __neq; ++__i) {
      _p[_slist1[__i]] += _p[_dlist1[__i]];
    }
  } else {
 B = B + wf ;
     }
   if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for general derivimplicit and KINETIC case */
    int __i, __neq = 4;
    double __state = C;
    double __primary_delta = (C + wf) - __state;
    double __dtsav = dt;
    for (__i = 0; __i < __neq; ++__i) {
      _p[_dlist1[__i]] = 0.0;
    }
    _p[_dlist1[2]] = __primary_delta;
    dt *= 0.5;
    v = NODEV(_pnt->node);
    _ode_matsol_instance1(_threadargs_);
    dt = __dtsav;
    for (__i = 0; __i < __neq; ++__i) {
      _p[_slist1[__i]] += _p[_dlist1[__i]];
    }
  } else {
 C = C + wf ;
     }
 _args[1] = _args[1] * d1 ;
   } }
 
static void _net_init(Point_process* _pnt, double* _args, double _lflag) {
    _args[1] = 1.0 ;
   _args[2] = t ;
   }
 
double Mgblock (  double _lv ) {
   double _lMgblock;
 _lMgblock = 1.0 / ( 1.0 + ( Mg / K0 ) * exp ( ( 0.001 ) * ( - z ) * delta * F * _lv / R / ( T + celsius ) ) ) ;
   
return _lMgblock;
 }
 
static double _hoc_Mgblock(void* _vptr) {
 double _r;
    _hoc_setdata(_vptr);
 _r =  Mgblock (  *getarg(1) );
 return(_r);
}
 
static int  rates (  double _lv ) {
   inf = ( _lv - gVDv0 ) * gVDst * gVI ;
   tau2 = ( tau2_0 + a2 * ( 1.0 - exp ( - b2 * _lv ) ) ) * q10_tau2 ;
   tau3 = ( tau3_0 + a3 * ( 1.0 - exp ( - b3 * _lv ) ) ) * q10_tau3 ;
   if ( tau1 / tau2 > .9999 ) {
     tau1 = .9999 * tau2 ;
     }
   if ( tau2 / tau3 > .9999 ) {
     tau2 = .9999 * tau3 ;
     }
    return 0; }
 
static double _hoc_rates(void* _vptr) {
 double _r;
    _hoc_setdata(_vptr);
 _r = 1.;
 rates (  *getarg(1) );
 return(_r);
}
 
static int _ode_count(int _type){ return 4;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
     _ode_spec1 ();
 }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 4; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 ();
 }
 
static void _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
 _ode_matsol_instance1(_threadargs_);
 }}

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  A = A0;
  B = B0;
  C = C0;
  gVD = gVD0;
 {
   Mgblock ( _threadargscomma_ v ) ;
   tau1 = tau1 * pow( Q10_tau1 , ( ( T0_tau - celsius ) / 10.0 ) ) ;
   q10_tau2 = pow( Q10_tau2 , ( ( T0_tau - celsius ) / 10.0 ) ) ;
   q10_tau3 = pow( Q10_tau3 , ( ( T0_tau - celsius ) / 10.0 ) ) ;
   tau = tauV * pow( Q10 , ( ( T0 - celsius ) / 10.0 ) ) ;
   rates ( _threadargscomma_ v ) ;
   wtau3 = 1.0 - wtau2 ;
   factor = - exp ( - tp / tau1 ) + wtau2 * exp ( - tp / tau2 ) + wtau3 * exp ( - tp / tau3 ) ;
   factor = 1.0 / factor ;
   A = 0.0 ;
   B = 0.0 ;
   C = 0.0 ;
   gVD = 0.0 ;
   wf = 1.0 ;
   }
  _sav_indep = t; t = _save;

}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _tsav = -1e20;
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
 initmodel();
}}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   i = ( wtau3 * C + wtau2 * B - A ) * ( gVI + gVD ) * Mgblock ( _threadargscomma_ v ) * ( v - e ) ;
   }
 _current += i;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 _g = _nrn_current(_v + .001);
 	{ _rhs = _nrn_current(_v);
 	}
 _g = (_g - _rhs)/.001;
 _g *=  1.e2/(_nd_area);
 _rhs *= 1.e2/(_nd_area);
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
double _dtsav = dt;
if (secondorder) { dt *= 0.5; }
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
 { error =  runge(_ninits, 4, _slist1, _dlist1, _p, &t, dt, state, &_temp1);
 if(error){fprintf(stderr,"at line 147 in file SynExp5NMDA.mod:\n\n"); nrn_complain(_p); abort_run(error);}
    if (secondorder) {
    int _i;
    for (_i = 0; _i < 4; ++_i) {
      _p[_slist1[_i]] += dt*_p[_dlist1[_i]];
    }}
  state();
 }}}
 dt = _dtsav;
}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(A) - _p;  _dlist1[0] = &(DA) - _p;
 _slist1[1] = &(B) - _p;  _dlist1[1] = &(DB) - _p;
 _slist1[2] = &(C) - _p;  _dlist1[2] = &(DC) - _p;
 _slist1[3] = &(gVD) - _p;  _dlist1[3] = &(DgVD) - _p;
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "/Users/hadasmanor/Documents/GitHubProjects/OOP_SimulationGUI/mods/SynExp5NMDA.mod";
static const char* nmodl_file_text = 
  "TITLE Triple-exp model of NMDAR has (HH-type gating) (temp. sensitivity) (voltage-dependent time constants) (desensitization)\n"
  "\n"
  "COMMENT\n"
  "This is a Triple-exponential model of an NMDAR \n"
  "that has a slow voltage-dependent gating component in its conductance\n"
  "time constants are voltage-dependent and temperature sensitive\n"
  "\n"
  "Mg++ voltage dependency from Spruston95 -> Woodhull, 1973 \n"
  "\n"
  "Desensitization is introduced in this model. Actually, this model has 4 differential equations\n"
  "becasue desensitization is solved analitically. It can be reduced to 3 by solving its A state analitically.\n"
  "For more info read the original paper. \n"
  "\n"
  "Keivan Moradi 2012\n"
  "\n"
  "ENDCOMMENT\n"
  "\n"
  "NEURON {\n"
  "	POINT_PROCESS Exp5NMDA\n"
  "	NONSPECIFIC_CURRENT i\n"
  "	RANGE tau1, tau2_0, a2, b2, wtau2, tau3_0, a3, b3, tauV, e, i, gVI, gVDst, gVDv0, Mg, K0, delta, tp, wf, tau_D1, d1\n"
  "	GLOBAL inf, tau2, tau3\n"
  "	THREADSAFE\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "	(nA) = (nanoamp)\n"
  "	(mV) = (millivolt)\n"
  "	(uS) = (microsiemens)\n"
  "	(mM) = (milli/liter)\n"
  "	(S)  = (siemens)\n"
  "	(pS) = (picosiemens)\n"
  "	(um) = (micron)\n"
  "	(J)  = (joules)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  ": Parameters Control Neurotransmitter and Voltage-dependent gating of NMDAR\n"
  "	tau1 = 1.69		(ms)	<1e-9,1e9>	: Spruston95 CA1 dend [Mg=0 v=-80 celcius=18] be careful: Mg can change these values\n"
  ": parameters control exponential rise to a maximum of tau2\n"
  "	tau2_0 = 3.97	(ms)\n"
  "	a2 = 0.70		(ms)\n"
  "	b2 = 0.0243		(1/mV)\n"
  "	wtau2= 0.65		<1e-9,1> : Hestrin90\n"
  "	\n"
  ": parameters control exponential rise to a maximum of tau3\n"
  "	tau3_0 = 41.62	(ms)\n"
  "	a3 = 34.69		(ms)\n"
  "	b3 = 0.01		(1/mV)\n"
  "	: Hestrin90 CA1 soma  [Mg=1 v=-40 celcius=30-32] the decay of the NMDA component of the EPSC recorded at temperatures above 30 degC \n"
  "	: the fast phase of decay, which accounted for 65%-+12% of the decay, had a time constant of 23.5-+3.8 ms, \n"
  "	: whereas the slow component had a time constant of 123-+83 ms.\n"
  "	: wtau2= 0.78 Spruston95 CA1 dend [Mg=0 v=-80 celcius=18] percentage of contribution of tau2 in deactivation of NMDAR\n"
  "	Q10_tau1 = 2.2			: Hestrin90\n"
  "	Q10_tau2 = 3.68			: Hestrin90 -> 3.5-+0.9, Korinek10 -> NR1/2B -> 3.68\n"
  "	Q10_tau3 = 2.65			: Korinek10\n"
  "	T0_tau	 = 35	(degC)	: reference temperature \n"
  "	: Hestrin90 CA1 soma  [Mg=1 v=-40 celcius=31.5->25] The average Q10 for the rising phase was 2.2-+0.5, \n"
  "	: and that for the major fast decaying phase was 3.5-+0.9\n"
  "	tp = 30			(ms)	: time of the peack -> when C + B - A reaches the maximum value or simply when NMDA has the peack current\n"
  "							: tp should be recalculated when tau1 or tau2 or tau3 changes\n"
  ": Parameters control desensitization of the channel\n"
  "	: these values are from Fig.3 in Varela et al. 1997\n"
  "	: the (1) is needed for the range limits to be effective\n"
  "	d1 = 0.2 	  	(1)		< 0, 1 >     : fast depression\n"
  "	tau_D1 = 2500 	(ms)	< 1e-9, 1e9 >\n"
  ": Parameters Control voltage-dependent gating of NMDAR\n"
  "	tauV = 7		(ms)	<1e-9,1e9>	: Kim11 \n"
  "							: at 26 degC & [Mg]o = 1 mM, \n"
  "							: [Mg]o = 0 reduces value of this parameter\n"
  "							: Because TauV at room temperature (20) & [Mg]o = 1 mM is 9.12 Clarke08 & Kim11 \n"
  "							: and because Q10 at 26 degC is 1.52\n"
  "							: then tauV at 26 degC should be 7 \n"
  "	gVDst = 0.007	(1/mV)	: steepness of the gVD-V graph from Clarke08 -> 2 units / 285 mv\n"
  "	gVDv0 = -100	(mV)	: Membrane potential at which there is no voltage dependent current, from Clarke08 -> -90 or -100\n"
  "	gVI = 1			(uS)	: Maximum Conductance of Voltage Independent component, This value is used to calculate gVD\n"
  "	Q10 = 1.52				: Kim11\n"
  "	T0 = 26			(degC)	: reference temperature \n"
  "	celsius 		(degC)	: actual temperature for simulation, defined in Neuron\n"
  ": Parameters Control Mg block of NMDAR\n"
  "	Mg = 1			(mM)	: external magnesium concentration from Spruston95\n"
  "	K0 = 4.1		(mM)	: IC50 at 0 mV from Spruston95\n"
  "	delta = 0.8 	(1)		: the electrical distance of the Mg2+ binding site from the outside of the membrane from Spruston95\n"
  ": The Parameter Controls Ohm haw in NMDAR\n"
  "	e = -0.7		(mV)	: in CA1-CA3 region = -0.7 from Spruston95\n"
  "}\n"
  "\n"
  "CONSTANT {\n"
  "	T = 273.16	(degC)\n"
  "	F = 9.648e4	(coul)	: Faraday's constant (coulombs/mol)\n"
  "	R = 8.315	(J/degC): universal gas constant (joules/mol/K)\n"
  "	z = 2		(1)		: valency of Mg2+\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "	v		(mV)\n"
  "	dt		(ms)\n"
  "	i		(nA)\n"
  "	g		(uS)\n"
  "	factor\n"
  "	wf\n"
  "	q10_tau2\n"
  "	q10_tau3\n"
  "	inf		(uS)\n"
  "	tau		(ms)\n"
  "	tau2	(ms)\n"
  "	tau3	(ms)\n"
  "	wtau3\n"
  "}\n"
  "\n"
  "STATE {\n"
  "	A		: Gating in response to release of Glutamate\n"
  "	B		: Gating in response to release of Glutamate\n"
  "	C		: Gating in response to release of Glutamate\n"
  "	gVD (uS): Voltage dependent gating\n"
  "}\n"
  "\n"
  "INITIAL { \n"
  "	Mgblock(v)\n"
  "	: temperature-sensitivity of the of NMDARs\n"
  "	tau1 = tau1 * Q10_tau1^((T0_tau - celsius)/10(degC))\n"
  "	q10_tau2 = Q10_tau2^((T0_tau - celsius)/10(degC))\n"
  "	q10_tau3 = Q10_tau3^((T0_tau - celsius)/10(degC))\n"
  "	: temperature-sensitivity of the slow unblock of NMDARs\n"
  "	tau  = tauV * Q10^((T0 - celsius)/10(degC))\n"
  "	\n"
  "	rates(v)\n"
  "	wtau3 = 1 - wtau2\n"
  "	: if tau3 >> tau2 and wtau3 << wtau2 -> Maximum conductance is determined by tau1 and tau2\n"
  "	: tp = tau1*tau2*log(tau2/(wtau2*tau1))/(tau2 - tau1)\n"
  "	\n"
  "	factor = -exp(-tp/tau1) + wtau2*exp(-tp/tau2) + wtau3*exp(-tp/tau3)\n"
  "	factor = 1/factor\n"
  "\n"
  "	A = 0\n"
  "	B = 0\n"
  "	C = 0\n"
  "	gVD = 0\n"
  "	wf = 1\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE state METHOD runge : derivimplicit : \n"
  "	: we found acceptable results with \"runge\" integration method\n"
  "	: However, M. Hines encouraged us to use \"derivimplicit\" method instead - which is slightly slower than runge - \n"
  "	: to avoid probable unstability problems\n"
  "\n"
  "	i = (wtau3*C + wtau2*B - A)*(gVI + gVD)*Mgblock(v)*(v - e)\n"
  "}\n"
  "\n"
  "DERIVATIVE state {\n"
  "	rates(v)\n"
  "	A' = -A/tau1\n"
  "	B' = -B/tau2\n"
  "	C' = -C/tau3\n"
  "	: Voltage Dapaendent Gating of NMDA needs prior binding to Glutamate Kim11\n"
  "	gVD' = ((wtau3*C + wtau2*B)/wf)*(inf-gVD)/tau\n"
  "	: gVD' = (inf-gVD)/tau\n"
  "}\n"
  "\n"
  "NET_RECEIVE(weight, D1, tsyn (ms)) {\n"
  "	INITIAL {\n"
  "	: these are in NET_RECEIVE to be per-stream\n"
  "	: this header will appear once per stream\n"
  "		D1 = 1\n"
  "		tsyn = t\n"
  "	}\n"
  "\n"
  "	D1 = 1 - (1-D1)*exp(-(t - tsyn)/tau_D1)\n"
  "	tsyn = t\n"
  "\n"
  "	wf = weight*factor*D1\n"
  "	A = A + wf\n"
  "	B = B + wf\n"
  "	C = C + wf\n"
  "\n"
  "	D1 = D1 * d1\n"
  "}\n"
  "\n"
  "FUNCTION Mgblock(v(mV)) {\n"
  "	: from Spruston95\n"
  "	Mgblock = 1 / (1 + (Mg/K0)*exp((0.001)*(-z)*delta*F*v/R/(T+celsius)))\n"
  "}\n"
  "\n"
  "PROCEDURE rates(v (mV)) { \n"
  "	inf = (v - gVDv0) * gVDst * gVI\n"
  "	\n"
  "	tau2 = (tau2_0 + a2*(1-exp(-b2*v)))*q10_tau2\n"
  "	tau3 = (tau3_0 + a3*(1-exp(-b3*v)))*q10_tau3\n"
  "	if (tau1/tau2 > .9999) {\n"
  "		tau1 = .9999*tau2\n"
  "	}\n"
  "	if (tau2/tau3 > .9999) {\n"
  "		tau2 = .9999*tau3\n"
  "	}\n"
  "}\n"
  ;
#endif
