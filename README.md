# PSO-Clustering-Wireless-Sensor-Networks
Master's project, University of Boumerdès, carried out as a pair project with Riham Boutaib.

Subject: Energy-efficient cluster head selection in a wireless sensor network with Particle Swarm Optimization,
compared with LEACH under the first-order radio model. 


## Problem

In a clustered sensor network, cluster heads relay and aggregate the data of their members, so they spend far
more energy than other nodes. LEACH rotates this role at random. Here a PSO chooses the heads from the
distance to the members and the residual energy, to extend the network lifetime.

## Model

First-order radio model, packet of k bits over a distance d:

```
Etx(k, d) = k*Eelec + k*efs*d^2      if d < d0
          = k*Eelec + k*emp*d^4      if d >= d0
Erx(k)    = k*Eelec
EDA       = m*k*EDA_bit              aggregation of m signals
d0        = sqrt(efs / emp)          about 87.7 m
```

| Parameter | Value |
|---|---|
| Field, nodes, base station | 100 x 100 m, 100 nodes, base station at 50, 50 |
| Initial energy | 80 normal nodes at 0.5 J, 20 advanced nodes at 1.0 J |
| Packet | 4000 bits |
| Eelec, efs, emp, EDA | 50 nJ/bit, 10 pJ/bit/m^2, 0.0013 pJ/bit/m^4, 5 nJ/bit/signal |
| LEACH | p = 0.05, epochs of 20 rounds |
| PSO | 30 particles, 100 iterations, w = 0.729, c1 = c2 = 1.494, velocity clamping |
| Fitness weights | alpha : beta : gamma = 0.5 : 0.3 : 0.2 |
| Re-clustering | every 5 rounds |
| Runs | 50 independent random deployments |

The MAC layer is ideal, with no collision or retransmission. The network is single hop and static, and a node is dead when its energy reaches 0.

## Algorithm

```
v = w*v + c1*r1*(pBest - x) + c2*r2*(gBest - x)
x = x + v
Fitness = alpha*Davg + beta*(1/Eres(CH)) + gamma*sigma_E
```

A particle holds K = 5 continuous positions, decoded into K distinct alive nodes.

### Implementation choices

- Each fitness term is normalised so that the weights are meaningful: Davg by the field diagonal, 1/Eres by the mean initial energy, sigma_E by the mean residual energy.
- sigma_E is the standard deviation of the residual energy predicted after the round. The standard deviation of the current energy does not depend on the particle, so it could not guide the search.
- The number of cluster heads is K = p * N = 5.
- The PSO runs again as soon as a cluster head dies, in addition to every 5 rounds.
- If no cluster head exists in a round, every alive node transmits directly to the base station.
- LEACH uses the epoch based eligibility set G, reset every 1/p rounds.
- The simulation stops when all nodes are dead, or at 5000 rounds.

## Metrics

First node death, last node death, their difference, total energy consumed, residual energy and its variance per round,
and packets delivered to the base station per round.

## Usage

```matlab
pso_wsn()                                   % report parameters, 50 runs
pso_wsn('nRuns', 5)                         % fewer runs
pso_wsn('nRuns', 3, 'swarm', 15, 'iters', 20)   % quick check
```

Output: `results/summary.csv` and, when graphics are available, `results/alive_nodes.png` and `results/residual_energy.png`.

The code uses no toolbox. It was run with GNU Octave 8.4. One PSO call with 30 particles and 100 iterations
takes about 4 s in Octave, so the full 50 runs are long there. Use fewer runs for a first look.

## Files

```
pso_wsn.m    complete simulation: parameters, radio model, LEACH, PSO, experiment and summary
```

## License

MIT
