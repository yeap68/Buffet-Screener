# app.py — Buffett-style Screener (Streamlit UI, iPhone-friendly deploy)
import math, statistics, urllib.request, json as jsonlib
from typing import Dict, List, Tuple, Optional
import streamlit as st
import pandas as pd
import plotly.express as px

SEC_BASE = "https://data.sec.gov"
# 👉 Put your email so SEC lets you in politely
HEADERS = {"User-Agent": "BuffetValueAnalyst/1.1 (Yeap Teck Hooi yeapteckhooi@gmail.com)"}

import time, urllib.error

def fetch_json(url: str) -> dict:
    last_err = None
    for i in range(5):
        try:
            req = urllib.request.Request(url, headers=HEADERS)
            with urllib.request.urlopen(req, timeout=30) as r:
                return jsonlib.loads(r.read().decode("utf-8"))
        except urllib.error.HTTPError as e:
            if e.code == 404:
                raise RuntimeError(f"URL 404 Not Found: {url}")
            if e.code in (403, 429, 503):
                time.sleep(1.5 * (i + 1))
                last_err = e
                continue
            raise
        except Exception as e:
            time.sleep(1.0 * (i + 1))
            last_err = e
            continue
    raise last_err

@st.cache_data(show_spinner=False, ttl=60*60)
def load_ticker_map() -> Dict[str, dict]:
    # Try primary XBRL endpoint
    sources = [
        f"{SEC_BASE}/api/xbrl/company_tickers.json",
        "https://www.sec.gov/files/company_tickers.json",  # fallback
    ]
    last_err = None
    for src in sources:
        try:
            data = fetch_json(src)
            if isinstance(data, dict) and data:
                return {v["ticker"].upper(): v for v in data.values()}
        except Exception as e:
            last_err = e
            continue

    # Final fallback: minimal hardcoded map so the app remains usable
    st.warning("SEC ticker map unavailable — using a limited built‑in map. Try again later for full coverage.")
    return {
        "AAPL": {"ticker": "AAPL", "cik_str": 320193},
        "MSFT": {"ticker": "MSFT", "cik_str": 789019},
        "AMZN": {"ticker": "AMZN", "cik_str": 1018724},
        "GOOGL": {"ticker": "GOOGL", "cik_str": 1652044},
        "BRK.B": {"ticker": "BRK.B", "cik_str": 1067983},
        "JNJ": {"ticker": "JNJ", "cik_str": 200406},
        "PG": {"ticker": "PG", "cik_str": 80424},
    }

def pick_usd_facts(facts, concept):
    try: units = facts["facts"]["us-gaap"][concept]["units"]
    except KeyError: return []
    out=[]
    for unit,items in units.items():
        if unit in ("USD","USD/shares","shares","pure"): out+=items
    return out

def rollup_annual(items):
    out={}
    for it in items:
        if it.get("fp") not in ("FY",None): continue
        if it.get("form") not in ("10-K","20-F","40-F","10-K/A","20-F/A","40-F/A"): continue
        fy, val = it.get("fy"), it.get("val")
        if fy and isinstance(val,(int,float)): out[int(fy)]=float(val)
    return out

def latest_ttm(facts, concept):
    try: units = facts["facts"]["us-gaap"][concept]["units"]
    except KeyError: return None
    q=[]
    for unit,items in units.items():
        for it in items:
            if it.get("fp") in ("Q1","Q2","Q3","Q4") and it.get("form") in ("10-Q","10-Q/A"):
                v=it.get("val"); 
                if isinstance(v,(int,float)): q.append((it.get("end"),float(v)))
    if len(q)<4: return None
    q.sort(key=lambda x:x[0]); return sum(v for _,v in q[-4:])

KEY = {
 "Revenue":"RevenueFromContractWithCustomerExcludingAssessedTax",
 "RevenueAlt":"SalesRevenueNet","GrossProfit":"GrossProfit","COGS":"CostOfGoodsAndServicesSold",
 "OperatingIncome":"OperatingIncomeLoss","PretaxIncome":"IncomeLossFromContinuingOperationsBeforeIncomeTaxesExtraordinaryItemsNoncontrollingInterest",
 "NetIncome":"NetIncomeLoss","CFO":"NetCashProvidedByUsedInOperatingActivities","Capex":"PaymentsToAcquirePropertyPlantAndEquipment",
 "SharesBasic":"WeightedAverageNumberOfSharesOutstandingBasic","SharesDiluted":"WeightedAverageNumberOfDilutedSharesOutstanding",
 "InterestExpense":"InterestExpense","Cash":"CashAndCashEquivalentsAtCarryingValue",
 "DebtLT":"LongTermDebtNoncurrent","DebtST":"DebtCurrent","Equity":"StockholdersEquity","Assets":"Assets","Segments":"NumberOfReportableSegments"
}

def choose_10y_plus_ttm(annual, ttm): 
    years=sorted(annual.keys())[-10:]; return years,[annual[y] for y in years],ttm
def invested_capital(assets,equity,cash,dl,ds):
    if equity is None: return None
    ic=(dl or 0)+(ds or 0)+equity-(cash or 0); return ic if ic and ic>0 else None
def interest_cov(ebit,ix): 
    if ebit is None or not ix: return None
    return ebit/abs(ix)

def pull_bundle(facts):
    rev_items = pick_usd_facts(facts, KEY["Revenue"]) or pick_usd_facts(facts, KEY["RevenueAlt"])
    rev_annual = rollup_annual(rev_items)
    rev_ttm = latest_ttm(facts, KEY["Revenue"]) or latest_ttm(facts, KEY["RevenueAlt"])
    def A(name): return rollup_annual(pick_usd_facts(facts,name))
    def T(name): return latest_ttm(facts,name)

    gpA, cogsA, opA, ptxA, niA, cfoA, capexA = A(KEY["GrossProfit"]),A(KEY["COGS"]),A(KEY["OperatingIncome"]),A(KEY["PretaxIncome"]),A(KEY["NetIncome"]),A(KEY["CFO"]),A(KEY["Capex"])
    shB, shD = A(KEY["SharesBasic"]), A(KEY["SharesDiluted"])
    ixA, cashA, dlA, dsA, eqA, asA, segA = A(KEY["InterestExpense"]),A(KEY["Cash"]),A(KEY["DebtLT"]),A(KEY["DebtST"]),A(KEY["Equity"]),A(KEY["Assets"]),A(KEY["Segments"])

    years, rev, rev_ttm = choose_10y_plus_ttm(rev_annual, rev_ttm)
    def align(d): return [d.get(y) for y in years]
    return {
      "years":years, "revenue":rev, "revenue_ttm":rev_ttm,
      "gross_profit":align(gpA),"cogs":align(cogsA),"operating_income":align(opA),
      "pretax_income":align(ptxA),"net_income":align(niA),
      "cfo":align(cfoA),"capex":align(capexA),
      "shares_basic":align(shB),"shares_diluted":align(shD),
      "interest_expense":align(ixA),"cash":align(cashA),"debt_lt":align(dlA),
      "debt_st":align(dsA),"equity":align(eqA),"assets":align(asA),"segments":align(segA),
      "ttm":{"revenue":rev_ttm,"gross_profit":T(KEY["GrossProfit"]),
             "operating_income":T(KEY["OperatingIncome"]),"net_income":T(KEY["NetIncome"]),
             "cfo":T(KEY["CFO"]),"capex":T(KEY["Capex"]),"interest_expense":T(KEY["InterestExpense"])}
    }

def compute_metrics(b):
    y=b["years"]; n=len(y)
    rev, gp, op, ni, cfo, capex = b["revenue"], b["gross_profit"], b["operating_income"], b["net_income"], b["cfo"], b["capex"]
    eq, cash, dl, ds, ix, assets = b["equity"], b["cash"], b["debt_lt"], b["debt_st"], b["interest_expense"], b["assets"]

    gm=[pct(gp[i],rev[i]) if i<len(gp) else None for i in range(n)]
    opm=[pct(op[i],rev[i]) if i<len(op) else None for i in range(n)]
    fcf=[(cfo[i] or 0)-(capex[i] or 0) if i<len(cfo) and i<len(capex) else None for i in range(n)]
    fcf_m=[pct(fcf[i],rev[i]) for i in range(n)]

    # ROIC ~ NOPAT / (Equity + Debt - Cash), tax ~ 21% fallback
    roic=[]
    for i in range(n):
        ic=invested_capital(assets[i] if i<len(assets) else None, eq[i] if i<len(eq) else None, 
                            cash[i] if i<len(cash) else None, dl[i] if i<len(dl) else None, ds[i] if i<len(ds) else None)
        ebit=op[i] if i<len(op) else None
        nopat=ebit*0.79 if ebit is not None else None
        roic.append(pct(nopat, ic) if (nopat is not None and ic is not None) else None)

    nd=(dl[-1] or 0)+(ds[-1] or 0)-(cash[-1] or 0)
    cov=interest_cov(op[-1], ix[-1] if ix else None)
    nd_to_ebit=pct(nd, op[-1]) if (op and op[-1] not in (None,0)) else None

    fcf_over_ni=[(fcf[i]/ni[i]) if (ni[i] not in (None,0) and fcf[i] is not None) else None for i in range(n)]
    rev_cagr_10=cagr(rev[0], rev[-1], max(1,y[-1]-y[0])) if rev and rev[0] and rev[-1] else None
    fcf_cagr_5=cagr(fcf[-6], fcf[-1], 5) if n>=6 and fcf[-6] and fcf[-1] else None
    s=safe(fcf); fcf_vol=(statistics.pstdev(s)/abs(statistics.mean(s))) if (len(s)>=5 and statistics.mean(s)!=0) else None

    return {"gm":gm,"opm":opm,"fcf":fcf,"fcf_margin":fcf_m,"roic":roic,
            "net_debt_to_ebit":nd_to_ebit,"interest_coverage":cov,
            "fcf_over_ni":fcf_over_ni,"rev_cagr_10y":rev_cagr_10,"fcf_cagr_5y":fcf_cagr_5,"fcf_vol":fcf_vol}

def owner_earnings_dcf(fcf_series, shares_series, price, discount, terminal_g, mos):
    sfcf=safe(fcf_series)
    if not sfcf: return {"iv_ps":None,"iv_ps_mos":None,"pass":False}
    fcf0=sfcf[-1]; sh=next((s for s in reversed(shares_series or []) if s), None)
    if not sh or sh<=0: return {"iv_ps":None,"iv_ps_mos":None,"pass":False}
    def cap_growth(g5,g10): return max(min(g5 or 0.0, (g10 or 0.0)*0.75), 0.0)
    g5=cagr(fcf_series[-6], fcf_series[-1], 5) if len(fcf_series)>=6 and fcf_series[-6] and fcf_series[-1] else 0.0
    g10=cagr(fcf_series[0], fcf_series[-1], max(1,len(fcf_series)-1)) if fcf_series and fcf_series[0] and fcf_series[-1] else 0.0
    g_base=cap_growth(g5,g10); g_bear=max(g_base-0.04,0.0); g_bull=g_base+0.04
    def pv_from(f, g, r, tg):
        horizon=5; lvl=f; pv=0.0; flows=[]
        for t in range(1,horizon+1):
            lvl*=(1+g); pv+=lvl/((1+r)**t); flows.append(lvl)
        if r<=tg: return None
        term=flows[-1]*(1+tg)/(r-tg); pv+=term/((1+r)**horizon); return pv
    ivb=pv_from(fcf0,g_base,discount,terminal_g); 
    ivbear=pv_from(fcf0,g_bear,discount+0.01,max(terminal_g-0.005,0.0))
    ivbull=pv_from(fcf0,g_bull,max(discount-0.01,0.06),min(terminal_g+0.005,0.035))
    if not ivb: return {"iv_ps":None,"iv_ps_mos":None,"pass":False}
    ivps=ivb/sh; ivps_mos=ivps*(1-mos)
    pass_mos = (price is not None and ivps_mos is not None and price <= ivps_mos)
    return {"iv_ps":float(ivps),"iv_ps_bear":float(ivbear/sh) if ivbear else None,
            "iv_ps_bull":float(ivbull/sh) if ivbull else None,
            "iv_ps_mos":float(ivps_mos),"pass":bool(pass_mos)}

def yahoo_price(ticker)->Optional[float]:
    try:
        req=urllib.request.Request(f"https://query1.finance.yahoo.com/v7/finance/quote?symbols={ticker}", headers=HEADERS)
        with urllib.request.urlopen(req, timeout=20) as r:
            data=jsonlib.loads(r.read().decode("utf-8"))
            return float(data["quoteResponse"]["result"][0]["regularMarketPrice"])
    except Exception:
        return None

def scorecard(m,b,mos_ok,min_roic=0.15,min_cov=6.0,max_nd_ebit=3.0,min_cash_conv=0.85):
    roic_med = statistics.median(safe(m["roic"])) if safe(m["roic"]) else None
    fcfni5 = statistics.mean(safe(m["fcf_over_ni"][-5:])) if m["fcf_over_ni"] else None
    gm=safe(m["gm"]); gm_decl=False
    if len(gm)>=4:
        try:
            d1=gm[-1]-gm[-2]; d2=gm[-2]-gm[-3]; d3=gm[-3]-gm[-4]
            gm_decl=((d1+d2+d3)/3)<=-0.02
        except: gm_decl=False
    segs=[s for s in b["segments"] if s is not None]; complex_biz=(max(segs) if segs else 1)>3
    return {
      "L1": not complex_biz,
      "L2": (roic_med is not None and roic_med>=min_roic and not gm_decl),
      "L3": (fcfni5 is not None and fcfni5>=min_cash_conv),
      "L4": ((m["net_debt_to_ebit"] is None or m["net_debt_to_ebit"]<max_nd_ebit) and (m["interest_coverage"] is None or m["interest_coverage"]>min_cov)),
      "L5": True,"L6": bool(mos_ok),"L7": True
    }

# ---------- UI ----------
st.set_page_config(page_title="Buffett Value Screener", layout="wide")
st.title("🧮 Buffett‑Style Value Screener")
st.caption("SEC XBRL (10Y+TTM) • Owner‑Earnings DCF • Adjustable MOS")

with st.sidebar:
    tickers_text = st.text_area("Tickers (comma‑separated)", "AAPL,MSFT,ADBE")
    mos = st.slider("Margin of Safety", 0.0, 0.7, 0.30, 0.01)
    discount = st.slider("Discount rate", 0.06, 0.14, 0.10, 0.005)
    terminal = st.slider("Terminal growth", 0.00, 0.035, 0.025, 0.001)
    min_roic = st.slider("Moat: 10‑yr median ROIC ≥", 0.05, 0.30, 0.15, 0.01)
    min_cov = st.slider("Interest coverage >", 2.0, 20.0, 6.0, 0.5)
    max_nd_ebit = st.slider("Net debt / EBIT <", 0.5, 8.0, 3.0, 0.1)
    min_cash_conv = st.slider("5‑yr avg FCF/NI ≥", 0.3, 1.5, 0.85, 0.05)
    run = st.button("Run screen", type="primary")

st.info("US filers with XBRL. Treat this as a **prefilter** and audit figures in the actual filings.", icon="ℹ️")

def analyze_one(t):
    tmap = load_ticker_map()
    if t not in tmap: return {"ticker": t, "error": "Not an SEC XBRL filer"}
    cik = str(tmap[t]["cik_str"])
    facts = get_company_facts(cik)
    b = pull_bundle(facts)
    m = compute_metrics(b)
    price = yahoo_price(t)
    dcf = owner_earnings_dcf(m["fcf"], b["shares_diluted"] or b["shares_basic"], price, discount, terminal, mos)
    sc = scorecard(m,b,dcf["pass"],min_roic=min_roic,min_cov=min_cov,max_nd_ebit=max_nd_ebit,min_cash_conv=min_cash_conv)
    return {"ticker":t,"bundle":b,"metrics":m,"price":price,"dcf":dcf,"score":sc}

if run:
    tickers=[x.strip().upper() for x in tickers_text.split(",") if x.strip()]
    rows=[]; results=[]
    prog=st.progress(0.0)
    for i,t in enumerate(tickers, start=1):
        results.append(analyze_one(t))
        prog.progress(i/len(tickers))
    for r in results:
        if "error" in r:
            rows.append({"Ticker":r["ticker"],"Error":r["error"]}); continue
        m, d, sc = r["metrics"], r["dcf"], r["score"]
        roic_med = statistics.median(safe(m["roic"])) if safe(m["roic"]) else None
        fcfni5 = statistics.mean(safe(m["fcf_over_ni"][-5:])) if m["fcf_over_ni"] else None
        gm_last = (safe(m["gm"][-1:])[0] if safe(m["gm"][-1:]) else None)
        opm_last = (safe(m["opm"][-1:])[0] if safe(m["opm"][-1:]) else None)
        rows.append({
          "Ticker":r["ticker"],"Price":r["price"],
          "IV_PS_Base":d.get("iv_ps"),"IV_PS_Bear":d.get("iv_ps_bear"),"IV_PS_Bull":d.get("iv_ps_bull"),
          "IV_PS_MOS":d.get("iv_ps_mos"),"MOS_Buy":d.get("pass"),
          "ROIC_10Y_Median":roic_med,"FCF/NI_5Y_Avg":fcfni5,
          "ND/EBIT_Last":m.get("net_debt_to_ebit"),"InterestCov_Last":m.get("interest_coverage"),
          "RevCAGR10Y":m.get("rev_cagr_10y"),"FCFCAGR5Y":m.get("fcf_cagr_5y"),"FCFVol":m.get("fcf_vol"),
          "GM_last":gm_last,"OPM_last":opm_last,
          "Pass_L1":sc["L1"],"Pass_L2":sc["L2"],"Pass_L3":sc["L3"],"Pass_L4":sc["L4"],
          "Pass_L5":sc["L5"],"Pass_L6":sc["L6"],"Pass_L7":sc["L7"],"AllPass":all(sc.values())
        })
    df=pd.DataFrame(rows)
    st.subheader("Results")
    st.dataframe(df, use_container_width=True)
    st.download_button("Download CSV", df.to_csv(index=False).encode("utf-8"), "screener_output.csv", "text/csv")
    shortlist = df[df["AllPass"]==True]
    st.success(f"Shortlist: {', '.join(shortlist['Ticker']) if not shortlist.empty else '—'}")
    st.markdown("---")
    st.subheader("Per‑ticker details")
    for r in results:
        if "error" in r:
            with st.expander(f"{r['ticker']} — Error"): st.error(r["error"]); continue
        b, m, d = r["bundle"], r["metrics"], r["dcf"]
        with st.expander(f"{r['ticker']} — audit & charts"):
            audit = pd.DataFrame({"Year":b["years"],"Revenue":b["revenue"],"FCF":m["fcf"],"ROIC":m["roic"]})
            st.dataframe(audit, use_container_width=True)
            chart_df = pd.DataFrame({"Year":b["years"],"FCF":m["fcf"],"ROIC":m["roic"]})
            c1,c2=st.columns(2)
            with c1: st.plotly_chart(px.line(chart_df, x="Year", y="FCF", markers=True, title="Free Cash Flow"), use_container_width=True)
            with c2: st.plotly_chart(px.line(chart_df, x="Year", y="ROIC", markers=True, title="ROIC"), use_container_width=True)
            st.markdown(f"**DCF (Base/Bear/Bull):** {d.get('iv_ps')} / {d.get('iv_ps_bear')} / {d.get('iv_ps_bull')}  •  **MOS IV:** {d.get('iv_ps_mos')}")
            st.caption("ROIC≈NOPAT/(Equity+Debt−Cash); ND/EBIT used as conservative proxy for ND/EBITDA.")
