"""Independent reconstruction, not the historical source. Round-end debit semantics."""
import numpy as np
from pathlib import Path
import json,time,csv,argparse,concurrent.futures
ROOT=Path(__file__).resolve().parent
CFG=json.loads((ROOT/'config.json').read_text())
PROT=list(CFG['protocol_stream_ids'])

def streams(seed,protocol):
    p=CFG['protocol_stream_ids'][protocol]
    return [np.random.default_rng(np.random.SeedSequence([seed,i])) for i in (0,1,p*10+2,p*10+3)]

def tx(d,b,c):
    q=c['radio']; d=np.asarray(d); return b*(q['elec']+np.where(d<np.sqrt(q['fs']/q['mp']),q['fs']*d*d,q['mp']*d**4))

def topology(c,seed,protocol):
    a,b,r,o=streams(seed,protocol); n=c['nodes']; area=np.array(c['area_m']); xy=a.uniform(0,1,(n,2))*area
    if c.get('layout')=='two_cluster':
        centers=np.array([[.25,.25],[.75,.75]])*area; xy=np.clip(centers[np.arange(n)%2]+a.normal(0,.09,(n,2))*area,0,area)
    h=c.get('heterogeneity',0); e=c['initial_energy_j']*(1+b.uniform(-h,h,n)); return xy,e,r,o

def geometry(ch,alive,dist,ds,e,c,multi):
    # Batched CH candidates; members attach nearest CH; strict progress excludes cycles.
    m,k=ch.shape; n=len(e); d=dist[:,ch].transpose(1,0,2); assignment=d.argmin(axis=2); near=d.min(axis=2); mask=np.broadcast_to(alive,(m,n)).copy(); np.put_along_axis(mask,ch,False,axis=1)
    load=np.stack([((assignment==i)&mask).sum(axis=1) for i in range(k)],axis=1)
    sd=ds[ch]; hops=np.full((m,k),k,dtype=int); rd=sd.copy()
    if multi:
        dd=dist[ch[:,:,None],ch[:,None,:]]; feasible=sd[:,None,:]<sd[:,:,None]-1e-12
        cost=tx(dd,c['payload_bits'],c)*(1+c['initial_energy_j']/np.maximum(e[ch[:,None,:]],1e-12)); cost=np.where(feasible,cost,np.inf)
        direct=tx(sd,c['payload_bits'],c); best=cost.argmin(axis=2); val=cost.min(axis=2); sel=val<direct
        hops=np.where(sel,best,k); rd=np.where(sel,np.take_along_axis(dd,best[:,:,None],axis=2)[:,:,0],sd)
    return assignment,near,mask,load,hops,rd

def objective(ch,alive,dist,ds,e,initial,duty,c,multi):
    _,near,mask,load,hops,rd=geometry(ch,alive,dist,ds,e,c,multi); k=ch.shape[1]
    diag=max(np.linalg.norm(c['area_m']),ds.max(),1e-12); ratio=np.clip(e[ch]/initial[ch],0,1)
    relay=np.stack([(hops==i).sum(axis=1) for i in range(k)],axis=1)
    terms=np.stack([(near*mask).sum(axis=1)/np.maximum(mask.sum(axis=1),1)/diag,rd.mean(axis=1)/diag,1-ratio.mean(axis=1),1-ratio.min(axis=1),load.std(axis=1)/np.maximum(load.mean(axis=1),1),duty[ch].mean(axis=1)/max(duty.max(),1),relay.max(axis=1)/max(k-1,1)],axis=1)
    return terms@np.array(c['fitness_weights'])

def choose(alive,e,initial,dist,ds,duty,last,r,protocol,c,rr,oo):
    ids=np.flatnonzero(alive); k=min(len(ids),max(1,int(np.floor(c['ch_fraction']*len(ids)+.5)))); multi=protocol.startswith('MH') or protocol=='WAOA-inspired'; epoch=int(np.ceil(1/c['ch_fraction']))*c['recluster_rounds']; pool=ids
    if not multi:
        pool=ids[(last[ids]<((r-1)//epoch)*epoch+1)]
        if len(pool)<k: pool=ids
    else:
        # Reconstruction decision: top 2K by energy / mean alive distance; not canonical distributed HEED.
        local=dist[np.ix_(ids,ids)].sum(axis=1)/max(len(ids)-1,1); score=(e[ids]/initial[ids])/(1+local/max(np.linalg.norm(c['area_m']),1e-12)); pool=ids[np.argsort(-score,kind='stable')[:min(len(ids),2*k)]]
    if protocol in ('LEACH','MH-HEED'):
        return np.sort(rr.choice(pool,k,replace=False)),[]
    a=c['agents']; it=c['iterations']; dim=len(pool); x=oo.random((a,dim)); v=np.zeros_like(x); pb=x.copy(); pf=np.full(a,np.inf); best=x[0].copy(); bf=np.inf; luc=np.ones(a); radius=np.full(a,np.sqrt(dim)); hist=[]
    for j in range(it+1):
        ch=np.sort(pool[np.argsort(-x,axis=1,kind='stable')[:,:k]],axis=1); f=objective(ch,alive,dist,ds,e,initial,duty,c,multi); b=f.argmin()
        if f[b]<bf: bf=float(f[b]); best=x[b].copy()
        hist.append(bf); improved=f<pf; pb[improved]=x[improved]; pf[improved]=f[improved]
        if j==it: break
        if 'GSO' in protocol:
            g=c['gso']; luc=g['retention']*luc+g['gain']/(1+f); delta=x[None,:,:]-x[:,None,:]; dd=np.linalg.norm(delta,axis=2); old=x.copy()
            for i in range(a):
                eligible=np.flatnonzero((dd[i]<radius[i])&(luc>luc[i])&(dd[i]>0)); radius[i]=np.clip(radius[i]+g['radius_gain']*(g['target_neighbors']-len(eligible)),0,np.sqrt(dim))
                if len(eligible):
                    w=luc[eligible]-luc[i]; z=oo.choice(eligible,p=w/w.sum()); x[i]=old[i]+g['step']*delta[i,z]/dd[i,z]
        elif protocol=='WAOA-inspired':
            # Local whale-inspired proxy; energy-progress route, not published WAOA.
            acoef=2*(1-j/max(it,1)); U=oo.random((a,dim)); A=2*acoef*U-acoef; C=2*oo.random((a,dim)); anchor=np.where(np.abs(A)<1,best,x[oo.integers(a,size=a)]); x=anchor-A*np.abs(C*anchor-x)
        else:
            p=c['pso']; w=p['w_start']+(p['w_end']-p['w_start'])*j/max(it-1,1); v=w*v+p['c1']*oo.random(x.shape)*(pb-x)+p['c2']*oo.random(x.shape)*(best-x)
            if protocol=='APSO-inspired': v+=oo.normal(0,.03,x.shape)*(1-j/max(it,1))
            v=np.clip(v,-p['vmax'],p['vmax']); x+=v
        x=np.clip(x,0,1)
    return np.sort(pool[np.argsort(-best,kind='stable')[:k]]),hist

def costs(ch,alive,dist,ds,e,c,multi,setup,retained_parent=None):
    n=len(e); b=c['payload_bits']; cb=c['control_bits']; q=c['radio']; ar,near,mask,load,hops,rd=geometry(ch[None,:],alive,dist,ds,e,c,multi); ar,near,mask,load,hops,rd=[x[0] for x in (ar,near,mask,load,hops,rd)]; k=len(ch)
    if retained_parent is not None:
        # Retain next hops until a new CH election; member losses change only receiver load.
        pos={int(v):i for i,v in enumerate(ch)}; hops=np.array([pos.get(int(retained_parent[z]),k) for z in ch]); rd=np.array([ds[z] if h==k else dist[z,ch[h]] for z,h in zip(ch,hops)])
    parent=np.full(n,-2,dtype=np.int16); role=np.zeros(n,dtype=np.int8); role[alive]=1; role[ch]=2; parent[mask]=ch[ar[mask]]; parent[ch]=np.where(hops==k,-1,ch[np.minimum(hops,k-1)]); relay=ch[np.unique(hops[hops<k])]; role[relay]=3
    txb=np.zeros(n,dtype=np.int32); rxb=txb.copy(); dab=txb.copy(); te=np.zeros(n); ce=np.zeros(n)
    txb[mask]=b; te[mask]=tx(near[mask],b,c); txb[ch]=b; te[ch]=tx(rd,b,c); rxb[ch]=load*b; dab[ch]=(load+1)*b
    # Each CH aggregates own sample + member samples + each upstream CH payload into one fixed-size packet.
    for h in hops[hops<k]: rxb[ch[h]]+=b; dab[ch[h]]+=b
    if setup:
        # CH advertises once at max member distance; every alive member listens to all CHs.
        for j,z in enumerate(ch):
            members=np.flatnonzero(mask&(ar==j)); radius=float(near[members].max()) if len(members) else 0.; ce[z]+=tx(radius,cb,c)+len(members)*cb*q['elec']; ce[members]+=tx(near[members],cb,c)+k*cb*q['elec']; ce[z]+=tx(rd[j],cb,c)
            if hops[j]<k: ce[ch[hops[j]]]+=cb*q['elec']
    debit=te+rxb*q['elec']+dab*q['da']+ce
    return debit,parent,role,txb,rxb,dab,ce

def simulate(task):
    phase,scenario,seed,protocol=task; c={**CFG,**CFG['stress'].get(scenario,{})}; start=time.perf_counter(); xy,e,rr,oo=topology(c,seed,protocol); initial=e.copy(); n=len(e); dist=np.linalg.norm(xy[:,None]-xy[None,:],axis=2); ds=np.linalg.norm(xy-np.array(c['sink_m']),axis=1); death=np.zeros(n,dtype=int); duty=np.zeros(n,dtype=int); last=np.full(n,-100000); ch=np.array([],int); r=1; next_setup=1; blocks=[]; convergence=[]; multi=protocol.startswith('MH') or protocol=='WAOA-inspired'; retained=None
    # Lossless event-block compression: constant costs between alive-set/setup changes.
    while r<=c['max_rounds'] and np.any(e>0):
        alive=e>0; setup=r>=next_setup or len(ch)==0 or not np.all(alive[ch]);
        if setup:
            ch,hist=choose(alive,e,initial,dist,ds,duty,last,r,protocol,c,rr,oo); duty[ch]+=1; last[ch]=r; next_setup=r+c['recluster_rounds']; retained=None
            if hist: convergence.append((r,hist))
        debit,parent,role,txb,rxb,dab,ce=costs(ch,alive,dist,ds,e,c,multi,setup,retained); retained=parent.copy()
        regular=debit-ce; remain=e-debit; ratios=np.divide(np.maximum(remain,0),regular,out=np.full(n,np.inf),where=(regular>0)&alive); soon=int(min(np.min(np.ceil(ratios[alive]-1e-12))+1,c['max_rounds']-r+1,next_setup-r)); span=max(1,soon)
        if np.any((remain<=0)&alive): span=1
        before=e.copy(); e=np.maximum(0,e-debit-(span-1)*regular); newly=alive&(e<=1e-13); e[newly]=0; death[newly]=r+span-1
        blocks.append((r,span,before.copy(),regular.copy(),ce.copy(),parent.copy(),role.copy(),txb.copy(),rxb.copy(),dab.copy()))
        r+=span
    path=ROOT/'outputs'/'logs'/phase/scenario/protocol; path.mkdir(parents=True,exist_ok=True); file=path/f'{seed}.npz'; keys=['start_round','span','energy_before','regular_debit','setup_debit','parent','role','tx_bits','rx_bits','aggregation_bits']; vals={k:np.asarray([z[i] for z in blocks]) for i,k in enumerate(keys)}; np.savez_compressed(file,**vals,xy=xy,initial=initial,death=death,ch_duty=duty,convergence_round=np.array([z[0] for z in convergence]),convergence=np.array([z[1] for z in convergence]).reshape(-1,c['iterations']+1))
    observed=death[death>0]; complete=len(observed)==n; t=np.sort(observed); h=(n+1)//2; fnd=float(t[0]) if len(t) else None; hnd=float(t[h-1]) if len(t)>=h else None; lnd=float(t[-1]) if complete else None; mean=float(t.mean()) if complete else None; jain=float(t.sum()**2/(n*(t*t).sum())) if complete else None
    return dict(phase=phase,scenario=scenario,seed=seed,protocol=protocol,nodes=n,FND=fnd,HND=hnd,LND=lnd,mean_lifetime=mean,Jain=jain,censored_nodes=n-len(t),runtime_s=time.perf_counter()-start,log=str(file.relative_to(ROOT)))

def tasks():
    return [(phase,'base',s,p) for phase,seeds in [('core',CFG['core_seeds']),('diagnostic',CFG['diagnostic_seeds']),('development',CFG['development_seeds'])] for s in seeds for p in PROT]+[('stress',sc,s,p) for sc in CFG['stress'] for s in CFG['stress_seeds'] for p in CFG['stress_protocols']]

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--workers',type=int,default=4); ap.add_argument('--phase',choices=['all','core','stress','development','diagnostic'],default='all'); args=ap.parse_args(); ts=[t for t in tasks() if args.phase=='all' or t[0]==args.phase]; out=ROOT/'outputs'/f'outcomes_{args.phase}.csv'; out.parent.mkdir(exist_ok=True)
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as ex, out.open('w',newline='') as f:
        w=None
        for i,row in enumerate(ex.map(simulate,ts),1):
            if w is None: w=csv.DictWriter(f,fieldnames=list(row)); w.writeheader()
            w.writerow(row); f.flush()
            if i%20==0 or i==len(ts): print(f'{i}/{len(ts)}',flush=True)
if __name__=='__main__': main()
