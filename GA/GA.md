# GA — Genetic Algorithm Pipeline

Source: John H. Holland, "Genetic Algorithms," *Scientific American* 267(1), July 1992, pp. 66–73.

## 1. Representation (the genetic code)

- Candidate solutions are encoded as fixed-length bit strings ("chromosomes"), e.g. `1101...`.
- Each bit position is a **gene**; the string's bit pattern is its **genotype**, and the behavior/solution it encodes is its **phenotype**.
- A classifier system generalizes this to condition→action rules: a string of `1`/`0`/`#` (don't-care) symbols, where the condition half matches input features and the action half specifies a response.
- Requirement: a small change to the genotype must produce a meaningful change in the phenotype (unlike, e.g., mutating raw FORTRAN source).

## 2. Initialize population

- Generate a population of random bit strings (Holland's classifier-system runs used up to ~8,000 rules).

## 3. Evaluate fitness

- Score every string against a problem-specific fitness measure (profitability, game payoff, error rate, engineering constraint satisfaction, etc.).
- No structural assumptions about the problem are required — fitness is a black-box measurement of solution quality.

## 4. Select parents

- Higher-fitness strings are given a higher chance to reproduce; low-fitness strings are discarded at each generation so population size stays constant.
- This selection pressure is what drives the *implicit parallelism* property: because a string belongs to every schema (bit-pattern region) it partially matches, sampling shifts toward high-fitness regions of the whole search space, not just toward individual good strings.

## 5. Crossover (recombination)

- Pick two parent strings and a random crossover point.
- Swap the tails beyond that point to produce two offspring, each combining building blocks from both parents.
- Offspring replace the lowest-fitness strings in the population (parents are not necessarily removed).
- Closely-spaced bits that jointly encode a useful trait ("compact building blocks") survive crossover more reliably than widely separated ones.

## 6. Mutation

- Flip a small fraction of bits at random (Holland cites roughly 1 in 10,000).
- Mutation alone rarely improves solutions; its role is to maintain diversity and prevent the population from converging to a uniform state that can no longer evolve.

## 7. (Optional) Inversion

- Occasionally reorders genes so that bits that were far apart become adjacent, converting a "loose" building block into a compact one that is less likely to be broken up by future crossovers.

## 8. Iterate

- Repeat steps 3–7 (evaluate → select → crossover → mutate) for many generations.
- Track whether the population's average and best fitness are improving; stop when a satisfactory solution is found, fitness plateaus, or a generation/time budget is exhausted.

## 9. (Extension) Rule competition for classifier systems

When evolving a system of many interacting rules rather than one best string:

- Every rule whose condition matches the current input "bids" to act; competing rules are weighted by past success, and the strongest wins control of the output.
- Winning rules are strengthened if the outcome is good, weakened otherwise (credit assignment can propagate back through a chain of rules for delayed payoffs).
- This produces default hierarchies: broad, simple rules emerge first as fallbacks; more specific exception rules evolve later and only override defaults when their conditions are met.

## 10. Practical notes from the article's case studies

- GAs excel at exploring large, nonlinear, poorly-understood search spaces (jet turbine design: >10^387 points; gas-pipeline control; network design; Prisoner's-Dilemma strategy discovery).
- They are weaker at fine local optimization — once near a good solution, hybridize with a local/hill-climbing method to polish it (as in the jet-turbine workflow: expert system → GA → refine).
- Balances exploration vs. exploitation automatically: crossover keeps testing new combinations of building blocks even while fit strings are being exploited, avoiding premature convergence on a single strategy.

## End-to-end pipeline summary

```
encode problem as bit strings
        |
        v
initialize random population
        |
        v
   +--> evaluate fitness --------------------+
   |            |                            |
   |            v                            |
   |     select parents (fitness-weighted)   |
   |            |                            |
   |            v                            |
   |     crossover -> offspring              |
   |            |                            |
   |            v                            |
   |     mutate (+ optional inversion)       |
   |            |                            |
   |            v                            |
   |     replace weakest strings -----------+
   |
   +-- repeat until stop condition (fitness target / plateau / generation limit)
        |
        v
  return best string(s) found
```
