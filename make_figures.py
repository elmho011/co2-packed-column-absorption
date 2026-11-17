"""Make the README charts from results.json.  Usage: python make_figures.py [results.json] [out_dir]"""
import json,sys,matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
R=json.load(open(sys.argv[1] if len(sys.argv)>1 else 'results.json')); out=sys.argv[2] if len(sys.argv)>2 else 'figures'
B,O='#2a78d6','#eb6834'; INK,MUT,GRID='#1f1f1e','#6b6a64','#e6e5e0'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.edgecolor':MUT,'axes.labelcolor':INK,
 'xtick.color':MUT,'ytick.color':MUT,'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,
 'grid.color':GRID,'grid.linewidth':0.8,'axes.titlesize':13,'axes.titleweight':'bold','axes.titlelocation':'left'})
def save(f,n): f.tight_layout(); f.savefig(f'{out}/{n}',dpi=160,facecolor='white'); plt.close(f)

# 1 parity
P=[p for p in R['pressure'] if p.get('dP_pred_ft')]
f,a=plt.subplots(figsize=(7,5))
mx=max(p['dP_pred_ft'] for p in P)*1.15; my=max(p['dP_meas_ft']+p['sd_ft'] for p in P)*1.08; m=max(mx,my)
a.plot([0,m],[0,m],color=MUT,lw=1.2,ls='--'); a.text(mx*0.97,mx*0.97-0.13,'1:1',color=MUT,ha='right')
a.plot([0,m/2],[0,m],color=MUT,lw=1,ls=':'); a.text(my/2-0.02,my*0.96,'2×',color=MUT,ha='right')
a.errorbar([p['dP_pred_ft'] for p in P],[p['dP_meas_ft'] for p in P],yerr=[p['sd_ft'] for p in P],fmt='o',ms=7,
  color=B,ecolor=B,elinewidth=1,capsize=2,mec='white',mew=1.5)
a.set_xlim(0,mx);a.set_ylim(0,my)
a.set_xlabel('GPDC predicted ΔP (in H₂O/ft)');a.set_ylabel('Measured ΔP (in H₂O/ft)')
a.set_title('Measured pressure drop runs ~2× the GPDC prediction')
save(f,'pressure_drop_parity.png')

# 2 KOGa time series
T=[x for x in R['kog_ts'] if x.get('k_cor') is not None and x.get('k_raw') is not None and 60<=x['t']<=680]; t=[x['t'] for x in T]
f,a=plt.subplots(figsize=(8,4.5))
a.plot(t,[x['k_cor'] for x in T],color=B,lw=2); a.plot(t,[x['k_raw'] for x in T],color=O,lw=2)
a.text(t[-1]+8,T[-1]['k_cor'],'Sensor-corrected',color=INK,va='center')
a.text(t[-1]+8,T[-1]['k_raw'],'Uncorrected',color=INK,va='center')
k=R['kog']; a.set_xlim(t[0],t[-1]*1.28)
a.set_ylim(0,max(x['k_cor'] for x in T)*1.2)
a.set_xlabel('Time into NaOH run (s)');a.set_ylabel('K$_{OG}$a  (kmol/h·m³·atm)')
a.set_title(f"K$_{{OG}}$a holds steady near {k['corrected']['KOGa']:.0f} kmol/(h·m³·atm) once corrected")
save(f,'kog_a_over_time.png')

# 3 driving force
D=R['driving']
f,a=plt.subplots(figsize=(8,4.5))
for key,c,lab in (('naoh',B,'0.076 M NaOH'),('water',O,'Water (model)')):
    z=[p['z']*1.219 for p in D[key]]; d=[p['drive'] for p in D[key]]
    a.plot(z,d,color=c,lw=2); a.text(z[-1]+0.02,d[-1],lab,color=INK,va='center')
a.set_xlim(0,1.219*1.28); a.set_ylim(0,None)
a.set_xlabel('Height up the packed bed (m)'); a.set_ylabel('Driving force, y − y*  (mol % CO₂)')
a.set_title('NaOH keeps the driving force high; water loses it')
save(f,'driving_force.png')
