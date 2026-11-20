# Ringkasan hasil (simulator rekonstruksi, 20 seed berpasangan, 50 node)

Sumber: outputs/*.csv. Selisih berpasangan dengan CI 95% (t) dan p Wilcoxon (asymptotic), Holm di dalam tiap tabel.

## A. Hasil utama (Table 3 baru)
| Protokol | FND | HND | LND | Jain | Mean lifetime |
|---|---|---|---|---|---|
| LEACH | 666.2 | 947.4 | 1336.7 | 0.9726 | 966.7 |
| MH-HEED | 901.0 | 1077.3 | 1121.8 | 0.9976 | 1061.5 |
| LEACH-GSO | 724.6 | 970.2 | 1221.6 | 0.9840 | 985.2 |
| LEACH-PSO | 710.0 | 981.4 | 1219.0 | 0.9855 | 986.6 |
| MH-HEED-GSO | 968.7 | 1091.7 | 1122.8 | 0.9989 | 1080.3 |
| MH-HEED-PSO | 975.9 | 1091.8 | 1123.0 | 0.9989 | 1080.5 |
| APSO-inspired | 717.9 | 978.5 | 1196.3 | 0.9863 | 986.1 |
| WAOA-inspired | 979.1 | 1091.8 | 1125.7 | 0.9989 | 1081.5 |

MH-HEED-PSO vs MH-HEED: FND +74.9 (8.3%), HND +14.55 (1.35%), mean lifetime +19.0 (1.79%), LND +1.2 (p=0.72, tidak berbeda).
MH-HEED-PSO vs MH-HEED-GSO: HND +0.15 (p=0.46), mean lifetime +0.2 (p=0.33).
Keluarga LEACH: PSO/GSO menaikkan HND (+34/+23) tetapi LND turun (~-117/-115): trade-off hanya di keluarga LEACH.

## B. Kontrol pada pool kandidat dan routing yang sama (E1, `E1_controls_*.csv`)
Selisih terhadap MH-HEED (pilih acak dari pool): HND / mean lifetime
- GREEDY (top-K skor HEED-like): -0.05 / +1.0 (tidak berbeda)
- Random search 88 evaluasi (RS): +13.95 / +20.1
- Hill-climbing 88 evaluasi (HC): +15.1 / +20.3
- Exhaustive (minimum eksak fungsi tujuan di pool): +12.55 / +19.3
- GSO: +14.4 / +18.8; PSO: +14.55 / +19.0
PSO vs RS: HND +0.6 (CI -1.7..2.9), mean lifetime -1.1 (CI -2.8..0.7). PSO vs EXH: HND +2.0 (CI -0.4..4.4).
Artinya manfaat berasal dari fungsi tujuan (pemilihan CH yang meminimalkan fungsi itu), bukan dari PSO/GSO: random search dengan anggaran sama setara, dan solusi eksak tidak lebih baik.

## C. Ablasi fungsi tujuan, MH-HEED-PSO (E2, `E2_ablation.csv`)
Hilangkan satu suku (bobot dinormalisasi ulang), selisih terhadap fungsi penuh:
- D_intra: HND -9.7, mean -8.1 (Holm p<0.001)
- CV_load: FND -46.4 (Holm 0.004), mean -4.9 (Holm 0.001)
- D_route: HND -4.1 (Holm 0.011), mean -1.2 (Holm 0.047)
- Bobot sama rata: HND -6.7 (Holm 0.003), mean -6.3 (Holm 0.001)
- P_meanE, P_minE, P_repeat, P_relay: tidak ada efek yang terdeteksi (Holm p>0.05)

## D. Pengulangan optimizer, 10 draw per topologi (E3)
SD antar draw dalam satu topologi vs SD antar topologi (untuk PSO, GSO, RS): HND ~4.0-4.1 vs ~5.7-6.0; mean lifetime ~2.3-2.4 vs ~5.5-6.1; FND ~20-21 vs ~10.5-12.6 (variasi optimizer mendominasi FND). MH-HEED acak: HND 5.9 vs 6.0, mean lifetime 3.5 vs 6.2, FND 30.1 vs 14.0.
Rata-rata 10 draw: PSO-GSO mean lifetime +0.77 (CI 0.33..1.20, Holm 0.004), kecil secara praktis; PSO-RS mean lifetime -0.61 (RS sedikit lebih baik).

## E. Metrik layanan (E4, E5)
Cakupan area (r_s = 15 m), AUC 0-1000 ronde: MH-HEED 0.9462, MH-HEED-PSO 0.9467, MH-HEED-GSO 0.9466, WAOA 0.9467, LEACH 0.9079.
PSO-MH-HEED: +0.0004 (CI 0.0001..0.0007); umur cakupan 90%: +17 ronde (Holm 0.019); PSO-LEACH: +0.039 dan +221.5 ronde. Efek sangat kecil dan bergantung pada r_s.
Beban terima per node (CV): LEACH 0.307, MH-HEED 0.116, PSO 0.087, GSO 0.089, WAOA 0.083. CV durasi CH: MH-HEED 0.311, PSO 0.152.
Catatan: tanpa batas jangkauan radio dan dengan re-elect CH langsung, setiap node hidup selalu punya rute ke sink; "sink-reachable fraction" = fraksi node hidup.

## F. Stress (Table 5 baru, HND rata-rata, 10 seed)
| Skenario | LEACH | MH-HEED | MH-HEED-PSO | WAOA |
|---|---|---|---|---|
| 25 node, 2k bit | 1806.9 | 2003.7 | 2032.7 | 2035.2 |
| 100 node, 6k bit | 660.2 | 759.9 | 759.3 | 759.2 |
| +-10% energi awal | 944.8 | 1079.5 | 1093.1 | 1094.4 |
| Sink tengah | 1137.0 | 1143.7 | 1159.4 | 1160.6 |
| Dua klaster | 993.2 | 1136.3 | 1140.2 | 1139.7 |
Sink tengah: PSO-LEACH +22.4 (10/10 seed), PSO-MH-HEED +15.7 (10/10), PSO-WAOA -1.2 (4/10). Skenario 100 node tanpa gain PSO. SD per sel ada di `Table5_sd.csv`.

## G. Runtime (E6, kondisi sama, satu proses, tanpa tulis log)
Rasio terhadap MH-HEED: GREEDY 0.75x, RS 1.8x, EXH 2.5x, HC 3.2x, PSO 5.6x, GSO 8.7x. Ini waktu simulasi di host, bukan biaya pada node.

## H. Belum ada / batasan
Pembanding faithful metode 2024+ (WAOA/APSO masih proxy), biaya trafik kontrol terpusat, validasi perangkat keras, sensitivitas independen untuk posisi sink, trafik, interval re-clustering dan bobot (hanya 5 skenario stress).
