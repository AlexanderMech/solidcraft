//! Butt and groove welds: one weld across a joint, filling the groove that is already modelled
//! in the parts.
//!
//! Not upstream: SolidCraft's own (#15132, #18500). The weld is built in two steps: the groove's
//! section, read from the parts ([`GrooveSection`]), and [`fill`], which fills it with the
//! user's options and returns the weld ([`GrooveWeld`]). Nothing here sees or changes a part.
//! For now the section is plain planar geometry the caller passes in (the two walls in the
//! plane square to the joint). Reading it from the parts, proved constant along a straight
//! joint, and the body swept along it are the next PR of the series (#18500).
//!
//! | Rust              | What it is                                                           |
//! |-------------------|----------------------------------------------------------------------|
//! | [`GrooveSection`] | The groove across the joint: its two walls and the plates' normal    |
//! | [`FillOptions`]   | The contour of each face and the penetration depth `s`               |
//! | [`Contour`]       | Flat, convex (sharp or blended toes) or concave                      |
//! | [`fill`]          | The weld: the groove filled, closed by its faces                     |
//! | [`GrooveWeld`]    | The weld's section, the groove's and the filled depth, its notes     |
//! | [`GrooveNote`]    | What the weld tool tells the user about a weld it built              |
//! | [`WeldFailure`]   | Why a weld is refused (shared by every weld type, `error::weld`)     |
//!
//! The options have no defaults here: the caller (ops) sets them (a flush face, a blend
//! overlap equal to the crest height).
//!
//! # The geometry
//!
//! The weld is generic: it fills whatever groove the walls bound, whatever its shape (V,
//! bevel, Y, U, J, the double-sided X and K, other angles, steps, rounds). The weld type is a
//! name for the symbol, not a recipe.
//!
//! Each wall is the groove face of one part as it cuts the section plane, a chain of lines and
//! arcs. It runs from its *toe* on the face side (where the groove face meets the part's outer
//! face) away from the face. Seen from the face side, `left` is on the left. A wall ends where
//! the groove ends: on the root side's outer face when the groove runs through the plate, or
//! where the two parts touch (both walls then end at the same point, the groove is closed).
//!
//! The *face* closes the weld between the two toes. Where the groove runs through the plate
//! the other side is closed between the two wall ends: flat on a weld from one side (its root,
//! always flush), or by the *second face* of a weld from both sides (a weld symbol above and
//! below the reference line: X, K, double-U, double-J, or any groove the user welds from both
//! sides), which follows the user's contour for that side. Each face is [`Contour::Flat`] (a
//! straight chord), [`Contour::Convex`] (an arc through the toes with its crest the given
//! height off the chord), [`Contour::ConvexBlended`] (the same crest, run out tangentially
//! past the toes) or [`Contour::Concave`] (an arc dipping the given depth below the chord).
//!
//! **Blended toes** (ISO 2553, "toes shall be blended smoothly"): a convex arc of radius `R`
//! between the toes and at each toe a concave arc of radius `r`, tangent to the chord's line
//! at its outer end and to the convex arc at the toe. The crest height is `h`, the blend runs
//! `o` past each toe. With half the chord `a` and `x = a + o`: `R + r = (x² + h²) / 2h`,
//! `R = a (R + r) / x`. `o = 0` is the plain convex face. Every arc is a polyline arc, the
//! chord between its ends and its bulge (`prism::path_export::bulge_arc`); none is more than
//! half a circle, so a convex face higher than half the distance between its toes is refused.
//! The blend runs along the chord's line, so it is built only where the plates' outer faces on
//! that side do ([`GrooveSection::plate_normal`] given, square to the chord); elsewhere it is
//! refused ([`WeldFailure::BlendedToesOffPlane`]).
//!
//! **Penetration depth `s`.** Where the plates' outer faces are parallel (a butt joint,
//! [`GrooveSection::plate_normal`] given), `s` is measured square to the plates from the face
//! side's outer face of the lower plate (the one whose toe lies deeper), and the weld's bottom
//! runs parallel to the plates: plates of unequal thickness or offset get `s` in the plate
//! whose surface is lower, the other plate more by the step between them. Where they are not
//! parallel (a bevel groove in a T- or corner joint) it is measured square to the face's chord
//! from the chord. The weld fills only that deep and the rest of the groove stays empty; the
//! second face then does not apply. Without `s` the whole groove is filled. An `s` deeper than
//! a closed groove fills the groove and returns a note; so does an `s` past the depth where a
//! groove that runs through opens on the root side (the thinner plate's thickness where the
//! plates' root sides are flush): full penetration. Neither is a refusal. [`GrooveWeld`]
//! returns the depth it filled.
//!
//! A concave face is built as asked and returns how deep it dips: below the chord, and below
//! each plate's outer face where the plates' normal is given (between plates in one plane that
//! is below the plate's surface: ISO 6520-1 no. 511). Its limit is closed form: the depth `F`
//! where it first touches another side of the section (a wall, the bottom, the other face);
//! the deepest that builds is `F - δ`, δ the depth margin (the larger of the section's
//! conservative allowance and the kernel's zone tolerance), refused past it
//! ([`WeldFailure::ConcaveFaceCrossesWall`]). With nothing touched, the half circle builds and
//! deeper is refused ([`WeldFailure::ConcaveFaceTooDeep`]). A face no higher or deeper than the
//! zone tolerance is built flat ([`GrooveNote::BuiltFlat`]); one past that but shallower than
//! the smallest arc the profile builds is refused ([`WeldFailure::FaceTooShallow`]).
//!
//! # Not covered
//!
//! Walls of lines and arcs only (conics and Béziers are refused). The section across one
//! station: a groove whose section changes along the joint is the caller's to find.

use std::f64::consts::{PI, TAU};

use crate::error::{Failure, Fit, FitMeasure, Reason, Result, SolidError};
use crate::kernel::Kernel;
use crate::math::{Point2, Vec2};
use crate::prism::path_export::{bulge_arc, Vertex};
use crate::prism::profile::{
    check_regions, smallest_arc_sweep, Loop2, ProfileRefusal, Region2, Seg2, MAX_ARC_RADIUS,
};

pub use crate::error::weld::{WeldFailure, WeldSide};

const OP: &str = "weld::groove";

/// How one side of the weld is closed between its toes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Contour {
    /// A straight chord between the toes: flush with plates in one plane.
    Flat,
    /// An arc through the toes, its crest `height` off the chord, the toes sharp.
    Convex {
        /// The crest's height above the chord (> 0).
        height: f64,
    },
    /// The convex arc with its toes blended: run out tangentially `overlap` past each toe.
    ConvexBlended {
        /// The crest's height above the chord (> 0).
        height: f64,
        /// How far the blend runs past each toe (>= 0; 0 is the plain convex face).
        overlap: f64,
    },
    /// An arc through the toes dipping `depth` below the chord.
    Concave {
        /// How deep the arc dips below the chord (> 0).
        depth: f64,
    },
}

/// One groove across its joint, as the weld tool reads it from the parts: the section of the
/// groove faces in the plane square to the joint.
#[derive(Debug, Clone, PartialEq)]
pub struct GrooveSection {
    /// The left part's groove face in the section plane, from its face-side toe away from the
    /// face: lines and arcs, each starting where the one before ends.
    pub left: Vec<Seg2>,
    /// The right part's groove face, the same way.
    pub right: Vec<Seg2>,
    /// The normal of the plates' outer faces on the face side, pointing out of the plates,
    /// where those faces are parallel (a butt joint, the plates in one plane or offset);
    /// `None` where they are not (a bevel groove in a T- or corner joint). It sets what the
    /// penetration depth is measured from.
    pub plate_normal: Option<Vec2>,
}

/// How the groove is filled: the user's options, every one set (ops gives the defaults).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FillOptions {
    /// The face's contour.
    pub face: Contour,
    /// The second face's contour for a weld from both sides (a symbol above and below the
    /// reference line); `None` for a weld from one side, whose root is flat. Used where the
    /// groove runs through the plate and is filled whole.
    pub second_face: Option<Contour>,
    /// The penetration depth `s` (see the module's docs for what it is measured from); `None`
    /// fills the whole groove.
    pub depth: Option<f64>,
}

/// What the weld tool tells the user about a weld it built.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum GrooveNote {
    /// The penetration depth asked for is deeper than the groove, which is closed where the
    /// parts touch: the whole groove is filled, no deeper.
    DeeperThanGroove {
        /// The depth asked for.
        asked: f64,
        /// The groove's depth.
        groove: f64,
    },
    /// The penetration depth asked for is past where the groove opens on the root side: the
    /// groove runs through and is filled whole (full penetration).
    DeeperThanPlate {
        /// The depth asked for.
        asked: f64,
        /// The depth where the groove opens on the root side (the thinner plate's thickness
        /// where the plates' root sides are flush).
        plate: f64,
    },
    /// A face asked for no higher or deeper than the kernel's zone tolerance: built flat (no
    /// arc builds that shallow).
    BuiltFlat {
        /// Which side.
        side: WeldSide,
        /// The height or depth asked for.
        asked: f64,
    },
    /// A concave contour dips below the chord between its toes.
    BelowToes {
        /// Which side.
        side: WeldSide,
        /// How deep below the chord.
        depth: f64,
        /// How deep the face's lowest point lies below each plate's outer face on that side,
        /// left then right: square to the plates where their normal is given, else below the
        /// chord (both `depth`).
        below: [f64; 2],
        /// `depth` over the groove's depth (the plate's thickness where the groove runs
        /// through).
        ratio: f64,
    },
}

/// A groove weld across its joint.
#[derive(Debug, Clone, PartialEq)]
pub struct GrooveWeld {
    /// The weld's section: one counter-clockwise loop, no holes.
    pub region: Region2,
    /// How deep the groove reaches, measured as the penetration depth is.
    pub groove_depth: f64,
    /// How deep the weld fills it, measured the same way: the penetration depth where the
    /// weld ends inside the groove, else the groove's depth.
    pub filled_depth: f64,
    /// Whether the groove runs through the plate (the walls end apart).
    pub through: bool,
    /// What to tell the user.
    pub notes: Vec<GrooveNote>,
}

/// The weld that fills `section` with `options`: the groove between its walls filled, closed
/// by its faces.
///
/// # Errors
///
/// [`SolidError::InvalidInput`] for walls that are empty, not lines and arcs, not connected,
/// not finite, that do not run away from the face, that cross each other, a plates' normal
/// that is not finite or points into the plates, or a contour or depth that is not finite and
/// positive. A face this construction does not build is refused as a [`WeldFailure`]: a
/// convex face past half a circle ([`WeldFailure::ConvexFaceTooHigh`]), a concave face that
/// reaches the section's other sides ([`WeldFailure::ConcaveFaceCrossesWall`]) or dips past
/// half a circle ([`WeldFailure::ConcaveFaceTooDeep`]), each with its limit as its [`Fit`], in
/// closed form; blended toes where the plates' outer faces do not run along the chord
/// ([`WeldFailure::BlendedToesOffPlane`]).
pub fn fill(k: &Kernel, section: &GrooveSection, options: &FillOptions) -> Result<GrooveWeld> {
    check_input(section, options)?;
    build(k, section, options)
}

/// The second face's contour: flat on a weld from one side.
fn second(options: &FillOptions) -> Contour {
    options.second_face.unwrap_or(Contour::Flat)
}

/// The weld.
fn build(k: &Kernel, section: &GrooveSection, options: &FillOptions) -> Result<GrooveWeld> {
    let (left, right) = (&section.left, &section.right);
    let (Some(l0), Some(r0), Some(lb), Some(rb)) = (
        left.first().map(Seg2::start),
        right.first().map(Seg2::start),
        left.last().map(Seg2::end),
        right.last().map(Seg2::end),
    ) else {
        return Err(SolidError::invalid_input(OP, "an empty wall"));
    };
    let eps = 1e-9 * scale(left.iter().chain(right));
    let chord = r0 - l0;
    let width = chord.length();
    if width <= eps {
        return Err(SolidError::invalid_input(OP, "the two toes at one point"));
    }
    let e = chord * (1.0 / width);
    let chord_up = left_normal(e);
    let under_chord = |p: Point2| -(p - l0).dot(chord_up);
    if under_chord(lb) <= eps || under_chord(rb) <= eps {
        return Err(SolidError::invalid_input(
            OP,
            "a wall that does not run away from the face (left and right swapped?)",
        ));
    }
    // what the penetration depth is measured from: the lower plate's outer face, square to
    // the plates, where they are parallel; else the chord
    let plates = match section.plate_normal {
        Some(n) => {
            // scaled by its larger part first, so neither its length's overflow nor its
            // underflow matters
            let m = n.x.abs().max(n.y.abs());
            let n = Vec2::new(n.x / m, n.y / m);
            let up = n * (1.0 / n.length());
            if up.dot(chord_up) <= 0.0 {
                return Err(SolidError::invalid_input(
                    OP,
                    "a plates' normal that points into the plates",
                ));
            }
            Some(up)
        }
        None => None,
    };
    let (origin, up) = match plates {
        Some(up) => (if (l0 - r0).dot(up) <= 0.0 { l0 } else { r0 }, up),
        None => (l0, chord_up),
    };
    let below = |p: Point2| -(p - origin).dot(up);
    let groove_depth = left
        .iter()
        .chain(right)
        .flat_map(|s| deepest_points(s, -up))
        .map(below)
        .fold(0.0_f64, f64::max);
    if groove_depth <= eps {
        return Err(SolidError::invalid_input(
            OP,
            "a groove of no depth below the lower plate's face",
        ));
    }
    let through = lb.distance(rb) > eps;
    // where a groove that runs through opens on the root side: past it the weld is full
    // penetration
    let full = if through {
        below(lb).min(below(rb))
    } else {
        groove_depth
    };
    let mut notes = Vec::new();

    // the walls as far as the weld fills them, and the closure at their far ends
    let (mut left_wall, mut right_wall) = (left.clone(), right.clone());
    let mut bottom: Vec<Seg2> = Vec::new();
    let mut cut = false;
    let mut filled_depth = groove_depth;
    if let Some(s) = options.depth {
        if s < full - eps {
            let (Some(l), Some(r)) = (
                cut_chain(left, s, origin, up, eps),
                cut_chain(right, s, origin, up, eps),
            ) else {
                return Err(SolidError::invalid_input(
                    OP,
                    "a penetration depth that cuts one wall only",
                ));
            };
            let (Some(le), Some(re)) = (l.last().map(Seg2::end), r.last().map(Seg2::end)) else {
                return Err(SolidError::invalid_input(
                    OP,
                    "a penetration depth above a toe",
                ));
            };
            left_wall = l;
            right_wall = r;
            if re.distance(le) > eps {
                bottom.push(Seg2::Line { a: re, b: le });
            }
            cut = true;
            filled_depth = s;
        } else if s > full + eps {
            notes.push(if through {
                GrooveNote::DeeperThanPlate {
                    asked: s,
                    plate: full,
                }
            } else {
                GrooveNote::DeeperThanGroove {
                    asked: s,
                    groove: groove_depth,
                }
            });
        }
    }
    // the two faces: what each would be, its own limits first
    let second_used = through && !cut;
    let zone = k.tolerances.zone_tol_3d;
    let mut shallow_options = *options;
    let face_floor;
    (shallow_options.face, face_floor) =
        shallow(WeldSide::Face, l0, r0, options.face, zone, &mut notes)?;
    let mut second_floor = 0.0;
    if second_used {
        let (contour, floor) =
            shallow(WeldSide::Second, rb, lb, second(options), zone, &mut notes)?;
        shallow_options.second_face = Some(contour);
        second_floor = floor;
    }
    let options = &shallow_options;
    if second_used {
        check_face(
            WeldSide::Second,
            rb,
            lb,
            second(options),
            plates.map(|n| -n),
            eps,
        )?;
    }
    check_face(WeldSide::Face, l0, r0, options.face, plates, eps)?;
    let second_path = if second_used {
        closure(rb, lb, second(options))
    } else {
        None
    };
    let face_path = closure(l0, r0, options.face);
    // a concave face: the deepest that touches none of the section's other sides
    let walls: Vec<Seg2> = right_wall.iter().chain(&left_wall).copied().collect();
    let mut limits = Vec::new();
    if let Contour::Concave { depth } = options.face {
        let mut others = walls.clone();
        if second_used {
            others.extend(
                second_path
                    .clone()
                    .unwrap_or_else(|| vec![Seg2::Line { a: rb, b: lb }]),
            );
        } else {
            others.extend(bottom.iter().copied());
        }
        let (fit, contact) = concave_limit(l0, r0, &others, eps);
        limits.push((WeldSide::Face, depth, fit, contact, face_floor));
    }
    if let (true, Contour::Concave { depth }) = (second_used, second(options)) {
        let mut others = walls;
        others.extend(
            face_path
                .clone()
                .unwrap_or_else(|| vec![Seg2::Line { a: l0, b: r0 }]),
        );
        let (fit, contact) = concave_limit(rb, lb, &others, eps);
        limits.push((WeldSide::Second, depth, fit, contact, second_floor));
    }
    // delta, the depth margin a face keeps from what it touches: the larger of the kernel's
    // zone tolerance and twice the near-miss tolerance of the crossing check (100 eps, eps
    // being the conservative allowance, 1e-9 of the section's size), so a face at its fit
    // stays clear of that check on any scale. A concave
    // face builds from its arc floor up to delta short of the touch (with nothing touched, up
    // to the half circle, which builds). Outside that, its fit is the nearest size that
    // builds: the floor below it, the deepest past it, the zone tolerance (built flat) where
    // the deepest is below the floor. Every fit builds when retried.
    let delta = zone.max(2.0 * near_miss(eps));
    for &(side, depth, touch, contact, floor) in &limits {
        let (upper, slack) = if contact {
            (touch - delta, 0.0)
        } else {
            (touch, eps)
        };
        if depth >= floor && depth <= upper + slack {
            continue;
        }
        let why = if depth < floor {
            WeldFailure::FaceTooShallow { side, size: depth }
        } else if contact {
            WeldFailure::ConcaveFaceCrossesWall { side, depth }
        } else {
            WeldFailure::ConcaveFaceTooDeep { side, depth }
        };
        let fit = if upper < floor {
            zone
        } else if depth < floor {
            floor
        } else {
            upper
        };
        return Err(Failure::new(OP, Reason::Weld(why))
            .with_fit(Fit::new(fit, FitMeasure::Distance, None))
            .into_error());
    }
    if second_used {
        bottom = second_path.ok_or_else(|| unbuilt(WeldSide::Second, rb, lb, second(options)))?;
    }
    let face = face_path.ok_or_else(|| unbuilt(WeldSide::Face, l0, r0, options.face))?;

    // the loop, clockwise: face, right wall down, bottom, left wall up
    let mut segs = face.clone();
    segs.extend(right_wall.iter().copied());
    segs.extend(bottom.clone());
    segs.extend(left_wall.iter().rev().map(Seg2::reversed));
    // the input check: no two sides of the loop cross. A face's own pieces are not checked
    // against each other (they are built to meet); each face is checked against the walls,
    // the bottom and the other face.
    let mut group = vec![1; face.len()];
    group.extend(right_wall.iter().map(|_| 0));
    group.extend(bottom.iter().map(|_| if second_used { 2 } else { 0 }));
    group.extend(left_wall.iter().map(|_| 0));
    if crosses(&segs, &group, near_miss(eps)) {
        return Err(SolidError::invalid_input(
            OP,
            "a section that crosses itself (walls that cross)",
        ));
    }
    let region = Region2 {
        outer: Loop2 { segs }.reversed(),
        holes: Vec::new(),
    };
    if let Err(refusal) = check_regions(std::iter::once(&region), eps.max(1e-12), zone) {
        return Err(SolidError::invalid_input(
            OP,
            match refusal {
                ProfileRefusal::LineOrArc => {
                    "a line or arc of the section that does not build (too short or too small)"
                }
                ProfileRefusal::Budget | ProfileRefusal::FreeForm => {
                    "a section that is not one closed loop"
                }
            },
        ));
    }
    if region.signed_area()?.sign() != Some(true) {
        return Err(SolidError::invalid_input(
            OP,
            "a section of no area (left and right swapped?)",
        ));
    }
    // a concave face's dip: below the chord, and below each plate's outer face on its side
    for (side, contour, closing, toes, out) in [
        (WeldSide::Face, options.face, &face, [l0, r0], plates),
        (
            WeldSide::Second,
            second(options),
            &bottom,
            [lb, rb],
            plates.map(|n| -n),
        ),
    ] {
        let used = side == WeldSide::Face || (through && !cut);
        let (true, Contour::Concave { depth }) = (used, contour) else {
            continue;
        };
        let below_toes = match out {
            Some(n) => {
                let lowest = |toe: Point2| {
                    closing
                        .iter()
                        .flat_map(|s| deepest_points(s, -n))
                        .map(|q| (toe - q).dot(n))
                        .fold(0.0_f64, f64::max)
                };
                [lowest(toes[0]), lowest(toes[1])]
            }
            None => [depth, depth],
        };
        notes.push(GrooveNote::BelowToes {
            side,
            depth,
            below: below_toes,
            ratio: depth / groove_depth,
        });
    }
    Ok(GrooveWeld {
        region,
        groove_depth,
        filled_depth,
        through,
        notes,
    })
}

fn check_input(section: &GrooveSection, options: &FillOptions) -> Result<()> {
    let eps = 1e-9 * scale(section.left.iter().chain(&section.right));
    for wall in [&section.left, &section.right] {
        if wall.is_empty() {
            return Err(SolidError::invalid_input(OP, "an empty wall"));
        }
        for s in wall {
            if !s.is_line_or_arc() {
                return Err(SolidError::invalid_input(
                    OP,
                    "a wall segment that is not a line or an arc",
                ));
            }
            let ok = match *s {
                Seg2::Line { a, b } => finite(a) && finite(b),
                Seg2::Arc {
                    center,
                    radius,
                    start,
                    sweep,
                } => {
                    finite(center)
                        && radius.is_finite()
                        && radius > 0.0
                        && start.is_finite()
                        && sweep.is_finite()
                        && sweep.abs() < TAU
                }
                _ => false,
            };
            if !ok {
                return Err(SolidError::invalid_input(
                    OP,
                    "a wall segment that is not finite",
                ));
            }
        }
        if wall
            .iter()
            .zip(wall.iter().skip(1))
            .any(|(a, b)| a.end().distance(b.start()) > eps.max(1e-12) * 10.0)
        {
            return Err(SolidError::invalid_input(
                OP,
                "a wall whose segments do not join",
            ));
        }
    }
    let positive = |v: f64| v.is_finite() && v > 0.0;
    if section
        .plate_normal
        .is_some_and(|n| !(n.x.is_finite() && n.y.is_finite()) || (n.x == 0.0 && n.y == 0.0))
    {
        return Err(SolidError::invalid_input(
            OP,
            "a plates' normal that is not finite and non-zero",
        ));
    }
    for contour in [options.face, second(options)] {
        let ok = match contour {
            Contour::Flat => true,
            Contour::Convex { height } => positive(height),
            Contour::ConvexBlended { height, overlap } => {
                positive(height) && overlap.is_finite() && overlap >= 0.0
            }
            Contour::Concave { depth } => positive(depth),
        };
        if !ok {
            return Err(SolidError::invalid_input(
                OP,
                "a contour height or depth that is not finite and positive",
            ));
        }
    }
    if options.depth.is_some_and(|s| !positive(s)) {
        return Err(SolidError::invalid_input(
            OP,
            "a penetration depth that is not finite and positive",
        ));
    }
    Ok(())
}

fn finite(p: Point2) -> bool {
    p.x.is_finite() && p.y.is_finite()
}

/// The size of the walls: their largest coordinate (arcs with their reach), at least 1.
fn scale<'a>(segs: impl Iterator<Item = &'a Seg2>) -> f64 {
    segs.fold(1.0_f64, |m, s| {
        let reach = match *s {
            Seg2::Arc { center, radius, .. } => center.x.abs().max(center.y.abs()) + radius,
            _ => 0.0,
        };
        m.max(s.start().max_dimension())
            .max(s.end().max_dimension())
            .max(reach)
    })
}

/// `v` turned a quarter turn counter-clockwise: the left of a run along `v`.
fn left_normal(v: Vec2) -> Vec2 {
    Vec2::new(-v.y, v.x)
}

/// The points of `s` that may lie furthest along `dir`: its ends, and an arc's point facing
/// `dir` where the arc passes it.
fn deepest_points(s: &Seg2, dir: Vec2) -> Vec<Point2> {
    let mut out = vec![s.start(), s.end()];
    if let Seg2::Arc { center, radius, .. } = *s {
        let p = center + dir * radius;
        if on_arc(s, p, 1e-12) {
            out.push(p);
        }
    }
    out
}

/// The path from `p` to `q` closing one side of the weld, its outside on the left of the run;
/// `None` for a convex or concave arc past half a circle.
fn closure(p: Point2, q: Point2, contour: Contour) -> Option<Vec<Seg2>> {
    let chord = q - p;
    let c = chord.length();
    let e = chord * (1.0 / c);
    let n = left_normal(e);
    let m = (p + q) * 0.5;
    let a = 0.5 * c;
    match contour {
        Contour::Flat => Some(vec![Seg2::Line { a: p, b: q }]),
        Contour::Convex { height } => sagitta_arc(p, q, height).map(|arc| vec![arc]),
        Contour::Concave { depth } => sagitta_arc(p, q, -depth).map(|arc| vec![arc]),
        Contour::ConvexBlended {
            height: h,
            overlap: o,
        } => {
            if o <= 0.0 {
                return sagitta_arc(p, q, h).map(|arc| vec![arc]);
            }
            let x = a + o;
            let sum = (x * x + h * h) / (2.0 * h);
            let big = a * sum / x;
            let small = sum - big;
            let at = |u: f64, v: f64| m + e * u + n * v;
            let c_main = at(0.0, h - big);
            let (c_left, c_right) = (at(-x, small), at(x, small));
            let inflect = |c_blend: Point2| c_main + (c_blend - c_main) * (big / sum);
            let (i_left, i_right) = (inflect(c_left), inflect(c_right));
            let (t_left, t_right) = (at(-x, 0.0), at(x, 0.0));
            Some(vec![
                Seg2::Line { a: p, b: t_left },
                arc_about(c_left, t_left, i_left, true)?,
                arc_about(c_main, i_left, i_right, false)?,
                arc_about(c_right, i_right, t_right, true)?,
                Seg2::Line { a: t_right, b: q },
            ])
        }
    }
}

/// The arc from `p` to `q` whose middle stands `sagitta` off the chord, on the left of the run
/// when positive, on its right when negative: the polyline arc of bulge `-2 sagitta / chord`
/// (`prism::path_export::bulge_arc`; a positive bulge turns counter-clockwise, bulging right
/// of the run). `None` past half a circle (`|sagitta|` more than half the chord).
fn sagitta_arc(p: Point2, q: Point2, sagitta: f64) -> Option<Seg2> {
    let bulge = -2.0 * sagitta / (q - p).length();
    bulge_arc(Vertex { at: p, bulge }, q)
}

/// A face too shallow for an arc: up to the zone tolerance it is built flat and a note says so.
/// Past it, a face lower than the smallest arc the profile builds (`floor`) does not build. A
/// convex face is refused here with that floor; a concave one is sized with it later, with
/// what it touches. The floor is the height at which the face's smallest arc sweeps more than
/// `smallest_arc_sweep` and its largest radius is below `MAX_ARC_RADIUS`: with `x = a + o`
/// (`o = 0` but for blended toes), the main arc sweeps `4 atan(h / x)` and each blend arc half
/// that, the main arc's radius is `a (x² + h²) / 2hx` and a blend arc's `o (x² + h²) / 2hx`:
/// the largest uses `m = max(a, o)`. The floor uses the sweep bound for
/// an infinite radius, dropping `smallest_arc_sweep`'s `zone / R` term: the sweep exceeds that
/// term for any `a > zone / 2`. The `(1 + 1e-9)` keeps a conservative allowance above it.
fn shallow(
    side: WeldSide,
    p: Point2,
    q: Point2,
    contour: Contour,
    zone: f64,
    notes: &mut Vec<GrooveNote>,
) -> Result<(Contour, f64)> {
    let (size, overlap) = match contour {
        Contour::Convex { height } => (height, 0.0),
        Contour::ConvexBlended { height, overlap } => (height, overlap.max(0.0)),
        Contour::Concave { depth } => (depth, 0.0),
        Contour::Flat => return Ok((contour, 0.0)),
    };
    if size <= zone {
        notes.push(GrooveNote::BuiltFlat { side, asked: size });
        return Ok((Contour::Flat, 0.0));
    }
    let a = 0.5 * p.distance(q);
    let x = a + overlap;
    let m = a.max(overlap);
    let b = MAX_ARC_RADIUS;
    let turn = if overlap > 0.0 { 0.5 } else { 0.25 };
    let floor = (x * (turn * smallest_arc_sweep(f64::MAX, zone)).tan())
        .max(x * m / (b + (b * b - m * m).max(0.0).sqrt()))
        * (1.0 + 1e-9);
    if size < floor && !matches!(contour, Contour::Concave { .. }) {
        return Err(
            Failure::new(OP, Reason::Weld(WeldFailure::FaceTooShallow { side, size }))
                .with_fit(Fit::new(floor, FitMeasure::Distance, None))
                .into_error(),
        );
    }
    Ok((contour, floor))
}

/// A face's own limits, before it is built: a convex face no higher than half a circle allows
/// (its fit `a + o`, half the toe distance plus the blend's overlap: the main arc's half sweep
/// `2 atan(h / (a + o))` reaches a quarter turn there), blended toes only where the plates'
/// outer faces on that side run along the chord (`out`, their normal, square to it).
fn check_face(
    side: WeldSide,
    p: Point2,
    q: Point2,
    contour: Contour,
    out: Option<Vec2>,
    eps: f64,
) -> Result<()> {
    let a = 0.5 * p.distance(q);
    let (height, overlap) = match contour {
        Contour::Convex { height } => (height, 0.0),
        Contour::ConvexBlended { height, overlap } => {
            let along = out.is_some_and(|n| (q - p).dot(n).abs() <= eps);
            if overlap > 0.0 && !along {
                return Err(Failure::new(
                    OP,
                    Reason::Weld(WeldFailure::BlendedToesOffPlane { side }),
                )
                .into_error());
            }
            (height, overlap.max(0.0))
        }
        Contour::Flat | Contour::Concave { .. } => return Ok(()),
    };
    let fit = a + overlap;
    if height > fit + eps {
        return Err(Failure::new(
            OP,
            Reason::Weld(WeldFailure::ConvexFaceTooHigh { side, height }),
        )
        .with_fit(Fit::new(fit, FitMeasure::Distance, None))
        .into_error());
    }
    Ok(())
}

/// Why a face whose limits hold still has no polyline arcs (its arcs a rounding past half a
/// circle at the limit): a convex face is too high, a concave one too deep, with its fit.
fn unbuilt(side: WeldSide, p: Point2, q: Point2, contour: Contour) -> SolidError {
    let a = 0.5 * p.distance(q);
    let (height, overlap) = match contour {
        Contour::Convex { height } => (height, 0.0),
        Contour::ConvexBlended { height, overlap } => (height, overlap.max(0.0)),
        Contour::Concave { depth } => {
            return Failure::new(
                OP,
                Reason::Weld(WeldFailure::ConcaveFaceTooDeep { side, depth }),
            )
            .with_fit(Fit::new(a, FitMeasure::Distance, None))
            .into_error();
        }
        Contour::Flat => return SolidError::invalid_input(OP, "a face that does not build"),
    };
    Failure::new(
        OP,
        Reason::Weld(WeldFailure::ConvexFaceTooHigh { side, height }),
    )
    .with_fit(Fit::new(a + overlap, FitMeasure::Distance, None))
    .into_error()
}

/// The depth where a concave arc from `p` to `q` (dipping on the right of the run) first
/// touches one of `others`, and whether it does (`false`: none is touched up to the half
/// circle, half the chord, which is then the limit). Closed form: the arcs through `p`
/// and `q` are a pencil of circles, each deeper one below the one before, so the first that
/// touches a side does so at one of its ends or where it is tangent to it. The circle of dip
/// `s` has its centre `m + n (R - s)`, `R = (a² + s²) / 2s` (`m` the chord's middle, `a` half
/// of it, `n` its normal on the outside); tangency with a line or a circle is a quadratic in
/// `s`, a point's circle is direct.
fn concave_limit(p: Point2, q: Point2, others: &[Seg2], eps: f64) -> (f64, bool) {
    let chord = q - p;
    let a = 0.5 * chord.length();
    let e = chord * (0.5 / a);
    let n = left_normal(e);
    let m = (p + q) * 0.5;
    let circle = |s: f64| {
        let r = (a * a + s * s) / (2.0 * s);
        (m + n * (r - s), r)
    };
    // within the arcs' reach: below the chord, between the toes
    let inside = |x: Point2| (x - m).dot(n) <= eps && (x - m).dot(e).abs() <= a + eps;
    let tol = near_miss(eps);
    let (mut best, mut contact) = (a, false);
    let mut take = |s: f64| {
        if s > 0.0 && s <= best {
            best = s;
            contact = true;
        }
    };
    let at_toe = |x: Point2| x.distance(p).min(x.distance(q)) <= tol;
    for seg in others {
        let (s0, s1) = (seg.start(), seg.end());
        // a side leaving a toe: the arc leaves it 2 atan(s / a) below the chord, so it meets a
        // side leaving at the angle t at s = a tan(t / 2) (a double root of the tangencies
        // below, so taken directly)
        for (toe, other) in [(p, q), (q, p)] {
            let Some(dir) = leaving(seg, toe, tol) else {
                continue;
            };
            let toward = (other - toe) * (1.0 / (2.0 * a));
            let angle = (-dir.dot(n)).atan2(dir.dot(toward));
            if angle > 0.0 && angle < PI {
                take(a * (0.5 * angle).tan());
            }
        }
        for x in [s0, s1] {
            let v = x - m;
            let under = v.dot(n);
            if x.distance(p).min(x.distance(q)) <= eps || under >= -eps || v.dot(e).abs() >= a {
                continue;
            }
            let h = (v.dot(v) - a * a) / (2.0 * under);
            // R - h, without the cancellation of the difference when h is large
            let r = (a * a + h * h).sqrt();
            take(if h > 0.0 { a * a / (r + h) } else { r - h });
        }
        match *seg {
            Seg2::Line { a: l0, b: l1 } => {
                let d = l1 - l0;
                let len = d.length();
                if len <= eps {
                    continue;
                }
                let t = d * (1.0 / len);
                let ln = left_normal(t);
                let (g, k) = ((m - l0).dot(ln), n.dot(ln));
                for sigma in [1.0, -1.0] {
                    // g + k (R - s) = sigma R, times 2s
                    for s in roots(-k - sigma, 2.0 * g, (k - sigma) * a * a) {
                        if s <= eps || s > a + eps {
                            continue;
                        }
                        let (c, _) = circle(s);
                        let x = c - ln * (c - l0).dot(ln);
                        let u = (x - l0).dot(t);
                        if (-tol..=len + tol).contains(&u) && inside(x) && !at_toe(x) {
                            take(s);
                        }
                    }
                }
            }
            Seg2::Arc { center, radius, .. } => {
                let w = m - center;
                let wn = w.dot(n);
                for sigma in [1.0, -1.0] {
                    // |c - center| = R + sigma radius, squared, times 2s
                    for s in roots(
                        -wn - sigma * radius,
                        w.dot(w) - a * a - radius * radius,
                        (wn - sigma * radius) * a * a,
                    ) {
                        if s <= eps || s > a + eps {
                            continue;
                        }
                        let (c, r) = circle(s);
                        let span = c - center;
                        let dist = span.length();
                        if dist <= eps {
                            continue;
                        }
                        let u = span * (1.0 / dist);
                        let x = [center + u * radius, center - u * radius]
                            .into_iter()
                            .min_by(|x, y| {
                                (x.distance(c) - r)
                                    .abs()
                                    .total_cmp(&(y.distance(c) - r).abs())
                            })
                            .unwrap_or(center);
                        if on_seg(seg, x, tol) && inside(x) && !at_toe(x) {
                            take(s);
                        }
                    }
                }
            }
            _ => {}
        }
    }
    (best, contact)
}

/// The unit direction a side runs in from `toe`, where it starts or ends there.
fn leaving(seg: &Seg2, toe: Point2, tol: f64) -> Option<Vec2> {
    let (start, end) = (
        seg.start().distance(toe) <= tol,
        seg.end().distance(toe) <= tol,
    );
    if !(start || end) {
        return None;
    }
    let d = match *seg {
        Seg2::Line { a, b } => {
            if start {
                b - a
            } else {
                a - b
            }
        }
        Seg2::Arc {
            start: s, sweep, ..
        } => {
            let at = if start { s } else { s + sweep };
            let travel = Vec2::new(-at.sin(), at.cos()) * sweep.signum();
            if start {
                travel
            } else {
                -travel
            }
        }
        _ => return None,
    };
    let len = d.length();
    (len > 0.0).then(|| d * (1.0 / len))
}

/// The real roots of `a x² + b x + c` (one of `b x + c` when `a` is zero).
fn roots(a: f64, b: f64, c: f64) -> Vec<f64> {
    if a.abs() <= 1e-15 * (b.abs() + c.abs()) {
        return if b == 0.0 { Vec::new() } else { vec![-c / b] };
    }
    let disc = b * b - 4.0 * a * c;
    if disc < 0.0 {
        return Vec::new();
    }
    // the root without cancellation, and the other from the product
    let r = -0.5 * (b + b.signum() * disc.sqrt());
    if r == 0.0 {
        return vec![0.0];
    }
    vec![r / a, c / r]
}

/// The arc about `center` from `p` to `q`, counter-clockwise when `ccw`, as the polyline arc
/// of its bulge `tan(sweep / 4)`; `None` past half a circle.
fn arc_about(center: Point2, p: Point2, q: Point2, ccw: bool) -> Option<Seg2> {
    let (dp, dq) = (p - center, q - center);
    let turn = (dq.y.atan2(dq.x) - dp.y.atan2(dp.x)).rem_euclid(TAU);
    let sweep = if ccw { turn } else { turn - TAU };
    bulge_arc(
        Vertex {
            at: p,
            bulge: (sweep / 4.0).tan(),
        },
        q,
    )
}

/// `chain` up to where it first reaches `s` below the line through `origin` square to `up`,
/// ending there; `None` when it never does.
fn cut_chain(chain: &[Seg2], s: f64, origin: Point2, up: Vec2, eps: f64) -> Option<Vec<Seg2>> {
    let below = |p: Point2| -(p - origin).dot(up);
    let mut out = Vec::new();
    for seg in chain {
        let hit = match *seg {
            Seg2::Line { a, b } => {
                let (da, db) = (below(a) - s, below(b) - s);
                (da * db <= 0.0 && da != db).then(|| da / (da - db))
            }
            Seg2::Arc {
                center,
                radius,
                start,
                sweep,
            } => {
                // the circle meets the line where below(center + r u(t)) = s
                let d = s - below(center);
                // below(center + r u) = below(center) - r (u . up): u . up = -d / r
                let k = -d / radius;
                if k.abs() > 1.0 {
                    None
                } else {
                    let base = up.y.atan2(up.x);
                    let off = k.acos();
                    [base + off, base - off]
                        .into_iter()
                        .filter_map(|ang| {
                            let t = if sweep >= 0.0 {
                                (ang - start).rem_euclid(TAU)
                            } else {
                                (start - ang).rem_euclid(TAU)
                            };
                            (t <= sweep.abs() + eps / radius).then_some(t / sweep.abs())
                        })
                        .fold(None, |m: Option<f64>, t| Some(m.map_or(t, |m| m.min(t))))
                }
            }
            _ => None,
        };
        match hit {
            Some(t) if t <= eps => {
                return (!out.is_empty()).then_some(out);
            }
            Some(t) if t < 1.0 - eps => {
                out.push(part(seg, t));
                return Some(out);
            }
            Some(_) => {
                out.push(*seg);
                return Some(out);
            }
            None => out.push(*seg),
        }
    }
    None
}

/// The first fraction `t` of a line or arc.
fn part(seg: &Seg2, t: f64) -> Seg2 {
    match *seg {
        Seg2::Line { a, b } => Seg2::Line {
            a,
            b: a + (b - a) * t,
        },
        Seg2::Arc {
            center,
            radius,
            start,
            sweep,
        } => Seg2::Arc {
            center,
            radius,
            start,
            sweep: sweep * t,
        },
        other => other,
    }
}

/// The near-miss tolerance of the crossing check: two sides closer than this meet.
fn near_miss(eps: f64) -> f64 {
    eps.max(1e-12) * 100.0
}

/// Whether the closed loop `segs` crosses or touches itself anywhere but where neighbours
/// join (`tol` round each join is not looked at). Neighbouring segments of one nonzero
/// `group` (a face's pieces) are one side: they are not checked against each other, and the
/// side joins its neighbours at its own ends.
fn crosses(segs: &[Seg2], group: &[usize], tol: f64) -> bool {
    // each segment's side, and each side's two ends
    let mut side = Vec::with_capacity(segs.len());
    let mut ends: Vec<(Point2, Point2)> = Vec::new();
    for (i, s) in segs.iter().enumerate() {
        let g = group.get(i).copied().unwrap_or(0);
        let joined = i > 0 && g != 0 && group.get(i - 1) == Some(&g);
        match ends.last_mut() {
            Some(last) if joined => last.1 = s.end(),
            _ => ends.push((s.start(), s.end())),
        }
        side.push(ends.len() - 1);
    }
    let n = ends.len();
    for (i, a) in segs.iter().enumerate() {
        for (j, b) in segs.iter().enumerate().skip(i + 1) {
            let (Some(&u), Some(&v)) = (side.get(i), side.get(j)) else {
                continue;
            };
            let Some(&(start, end)) = ends.get(u) else {
                continue;
            };
            if u == v {
                continue;
            }
            let next = v == u + 1;
            let wrap = u == 0 && v + 1 == n;
            let joins: Vec<Point2> = match (next, wrap) {
                (true, true) => vec![end, start],
                (true, false) => vec![end],
                (false, true) => vec![start],
                (false, false) => Vec::new(),
            };
            if meets(a, b, tol)
                .into_iter()
                .any(|p| joins.iter().all(|&v| p.distance(v) > tol))
            {
                return true;
            }
        }
    }
    false
}

/// Where two lines or arcs meet (where they overlap, the ends of the overlap).
fn meets(a: &Seg2, b: &Seg2, tol: f64) -> Vec<Point2> {
    let candidates = match (*a, *b) {
        (Seg2::Line { a: p, b: q }, Seg2::Line { a: r, b: s }) => line_line(p, q, r, s, tol),
        (Seg2::Line { a: p, b: q }, Seg2::Arc { center, radius, .. })
        | (Seg2::Arc { center, radius, .. }, Seg2::Line { a: p, b: q }) => {
            line_circle(p, q, center, radius)
        }
        (
            Seg2::Arc {
                center: c1,
                radius: r1,
                ..
            },
            Seg2::Arc {
                center: c2,
                radius: r2,
                ..
            },
        ) => {
            if c1.distance(c2) <= tol && (r1 - r2).abs() <= tol {
                // one circle: the arcs meet where an end of one lies on the other
                vec![a.start(), a.end(), b.start(), b.end()]
            } else {
                circle_circle(c1, r1, c2, r2)
            }
        }
        _ => Vec::new(),
    };
    candidates
        .into_iter()
        .filter(|&p| on_seg(a, p, tol) && on_seg(b, p, tol))
        .collect()
}

/// Whether `p` lies on the line or arc `s`, within `tol`.
fn on_seg(s: &Seg2, p: Point2, tol: f64) -> bool {
    match *s {
        Seg2::Line { a, b } => {
            let d = b - a;
            let len2 = d.length_squared();
            if len2 == 0.0 {
                return p.distance(a) <= tol;
            }
            let t = (p - a).dot(d) / len2;
            let near = a + d * t.clamp(0.0, 1.0);
            p.distance(near) <= tol
        }
        Seg2::Arc { center, radius, .. } => {
            (p.distance(center) - radius).abs() <= tol && on_arc(s, p, tol)
        }
        _ => false,
    }
}

/// Whether the direction of `p` from an arc's centre lies within its sweep (`tol` along it).
fn on_arc(s: &Seg2, p: Point2, tol: f64) -> bool {
    let Seg2::Arc {
        center,
        radius,
        start,
        sweep,
    } = *s
    else {
        return false;
    };
    let d = p - center;
    let ang = d.y.atan2(d.x);
    let t = if sweep >= 0.0 {
        (ang - start).rem_euclid(TAU)
    } else {
        (start - ang).rem_euclid(TAU)
    };
    let slack = tol / radius;
    t <= sweep.abs() + slack || t >= TAU - slack
}

fn line_line(p: Point2, q: Point2, r: Point2, s: Point2, tol: f64) -> Vec<Point2> {
    let (d1, d2) = (q - p, s - r);
    let den = d1.x * d2.y - d1.y * d2.x;
    let scale = d1.length() * d2.length();
    if den.abs() <= 1e-14 * scale {
        // parallel: they meet only if they lie on one line, at the ends of the overlap
        let off = (r - p).x * d1.y - (r - p).y * d1.x;
        if off.abs() <= tol * d1.length() {
            return vec![p, q, r, s];
        }
        return Vec::new();
    }
    let w = r - p;
    let t = (w.x * d2.y - w.y * d2.x) / den;
    vec![p + d1 * t]
}

fn line_circle(p: Point2, q: Point2, c: Point2, r: f64) -> Vec<Point2> {
    let d = q - p;
    let f = p - c;
    let a = d.dot(d);
    if a == 0.0 {
        return Vec::new();
    }
    let b = 2.0 * f.dot(d);
    let cc = f.dot(f) - r * r;
    let disc = b * b - 4.0 * a * cc;
    if disc < 0.0 {
        // a near miss within rounding counts as touching at the foot of the centre
        let t = -b / (2.0 * a);
        return vec![p + d * t];
    }
    let root = disc.sqrt();
    vec![
        p + d * ((-b - root) / (2.0 * a)),
        p + d * ((-b + root) / (2.0 * a)),
    ]
}

fn circle_circle(c1: Point2, r1: f64, c2: Point2, r2: f64) -> Vec<Point2> {
    let d = c2 - c1;
    let dist = d.length();
    if dist == 0.0 {
        return Vec::new();
    }
    let a = (dist * dist + r1 * r1 - r2 * r2) / (2.0 * dist);
    let h2 = r1 * r1 - a * a;
    let u = d * (1.0 / dist);
    let mid = c1 + u * a;
    if h2 < 0.0 {
        // a near miss within rounding counts as touching on the line of centres
        return vec![mid];
    }
    let h = h2.sqrt();
    let n = left_normal(u);
    vec![mid + n * h, mid - n * h]
}
