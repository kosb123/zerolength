const force = [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  240000,  0,  0,  0,  0,  -60000,  0,  0,  0,  -180000,  0,  0,  0,  0,  0,  0,]

const localdistributedloads = [
  { member_id: 1, wy: -4000, wz: 0 },
  { member_id: 2, wy: -4000, wz: 0 }
];

const ed = new Array(12).fill(0);
const edp = new Array(12).fill(0);

function AddLoads(force){
    localdistributedloads.forEach(load => {


//----- loads due to uniformly distributed load on element
     for (n = 0; n < ne; n++) {
           ElemStiff(n)
           i1 = noc[n][0] - 1
           i2 = noc[n][1] - 1
           x21 = x[i2][0] - x[i1][0]
           y21 = x[i2][1] - x[i1][1]
           z21 = x[i2][2] - x[i1][2]
           el = Math.sqrt(x21 * x21 + y21 * y21 + z21 * z21)
           ed[0] = 0
           ed[1] = udl[n][0] * el / 2
           ed[2] = udl[n][1] * el / 2
           ed[3] = 0
           ed[4] = -udl[n][1] * el * el / 12
           ed[5] = udl[n][0] * el * el / 12
           ed[6] = 0
           ed[7] = ed[1]
           ed[8] = ed[2]
           ed[9] = 0
           ed[10] = -ed[4]
           ed[11] = -ed[5]
           for (i = 0; i < 12; i++) {
              edp[i] = 0
              for (k = 0; k < 12; k++) {
                 edp[i] = edp[i] + alambda[k][i] * ed[k]
              }
           }
           for (i = 0; i < 6; i++) {
              force[6 * i1 + i] = force[6 * i1 + i] + edp[i]
              force[6 * i2 + i] = force[6 * i2 + i] + edp[i + 6]
           }
        }
     }
