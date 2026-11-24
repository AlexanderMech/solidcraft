//! Groove welds (#18500): the joint of a groove (butt) weld, read from the two parts the user
//! picked: each part's groove faces as one chain from the toe ([`chain`]), the joint's direction
//! and frame ([`joint_frame`]), and the weld that fills the groove ([`groove_weld`]).
//!
//! SolidCraft's own: solidcore and solidcore-ops hold no weld (bherbruck on #18409). Every read
//! here is an ops call; this module adds the weld's rules and refusals. Nothing writes to a part.
//!
//! | Rust              | What it is                                                         |
//! |-------------------|--------------------------------------------------------------------|
//! | [`GroovePick`]    | One part's side of the joint: its groove faces and its outer face   |
//! | [`JointPart`]     | Which of the two parts something is about                           |
//! | [`Chain`]         | A part's groove faces in order from the toe, and its toe edges      |
//! | [`chain`]         | The chain of one pick                                               |
//! | [`joint_frame`]   | The joint's frame: a plane whose normal is the joint's direction    |
//! | [`groove_weld`]   | The weld: the groove's walls, its fill and its body                  |
//! | [`GrooveWeld`]    | What [`groove_weld`] gives                                          |
//! | [`GrooveRefusal`] | Why a joint is refused                                              |
//!
//! # Not covered
//!
//! T joints, whose toe on the through part is no edge, are refused ([`GrooveRefusal::TJoint`]);
//! their toe and the seam mark come in a later PR (brett19 on #18500). A joint whose section
//! changes along it is refused ([`GrooveRefusal::SectionChanges`]): its stretches of one section,
//! each built on its own, come next. So are parts that do not run alongside each other over the
//! same length.

use std::f64::consts::FRAC_PI_2;

use solidcore_ops::prism::{ChannelFill, FillFailure, FillOptions, PlanarPiece, extrude_planar_regions, fill_channel, planar_edges};
use solidcore_ops::{Body, Order, affine, sweep};
use solidcore_record::{Plane, Seg2, Vec2, Vec3};

/// Which of the two parts of a joint something is about.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum JointPart {
    /// The left part (the first one passed).
    Left,
    /// The right part.
    Right,
}

/// Why a groove joint is refused.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum GrooveRefusal {
    /// Input the joint cannot take: no groove face picked, or a face the part does not have.
    #[error("{0}")]
    Invalid(&'static str),
    /// The groove faces picked on a part do not form one chain of faces that share edges, from
    /// the toe to the root: record face `face` is off the chain (not joined to it, or a third
    /// neighbour of a face on it).
    #[error("the {part:?} part's groove faces are not one chain: face {face} is off it")]
    PicksNotOneChain { part: JointPart, face: usize },
    /// A part's outer face `face` (on the arrow side) does not meet its groove faces at one end
    /// of their chain: there is no toe.
    #[error("the {part:?} part's outer face {face} does not meet its groove faces at a toe")]
    ArrowFaceNotAtToe { part: JointPart, face: usize },
    /// A part's outer face is one of its groove faces: a T joint, whose toe on the through part
    /// is no edge. Not built yet.
    #[error("the {part:?} part's outer face is one of its groove faces (a T joint)")]
    TJoint { part: JointPart },
    /// A part's outer face `face` is not a plane.
    #[error("the {part:?} part's outer face {face} is not a plane")]
    OuterFaceNotPlane { part: JointPart, face: usize },
    /// Face `face` of a part (a groove face or its outer face) is not proved straight along the
    /// joint's direction within the zone tolerance ([`Body::face_straying`]).
    #[error("face {face} of the {part:?} part is not straight along the joint")]
    GrooveFaceNotStraight { part: JointPart, face: usize },
    /// Face `face` of a part is not read as the same lines and arcs at both of its ends in the
    /// joint's frame: the joint's section changes along it, or the face is no line or arc there.
    /// Only a joint of one section all along is built yet.
    #[error("face {face} of the {part:?} part does not keep one section along the joint")]
    SectionChanges { part: JointPart, face: usize },
    /// The two parts do not run alongside each other over the same length of the joint (their
    /// toes differ in length or are shifted along it). Not built yet.
    #[error("the parts do not run alongside each other over the same length of the joint")]
    PartsNotAlongside,
    /// The groove cannot be filled as asked ([`FillFailure`]); `fit` is the nearest size that
    /// builds, where there is one.
    #[error("the groove cannot be filled as asked: {why}{}", fit.map_or(String::new(), |f| format!(" (at most {f} builds)")))]
    Fill { why: FillFailure, fit: Option<f64> },
    /// An ops call failed.
    #[error(transparent)]
    Ops(#[from] solidcore_ops::Error),
}

type Result<T> = std::result::Result<T, GrooveRefusal>;

/// One part's side of a groove joint.
#[derive(Debug, Clone, Copy)]
pub struct GroovePick<'a> {
    /// The part (only read).
    pub part: &'a Body,
    /// Its record order: the faces below are record indices in it.
    pub order: &'a Order,
    /// Its groove faces (bevel, root face, J or U radius, square wall), in any order.
    pub faces: &'a [usize],
    /// Its plate face on the arrow side: the side the weld is made from.
    pub outer: usize,
}

/// A part's groove faces in order from the toe.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Chain {
    /// The groove faces (record indices), the one beside the outer face first, each next one
    /// sharing an edge with the one before.
    pub faces: Vec<usize>,
    /// The toe: the edges (record indices) the first face shares with the outer face.
    pub toe: Vec<usize>,
}

/// The picked groove faces of `pick` as one chain from the toe, read from the part's topology
/// only: which faces share an edge ([`Body::rim_edges`] of each face). The order of the picks
/// changes nothing. `part` names the part in a refusal.
///
/// # Errors
/// [`GrooveRefusal::Invalid`] for no groove face, or a face past the part's faces;
/// [`GrooveRefusal::TJoint`] when the outer face is among the groove faces;
/// [`GrooveRefusal::PicksNotOneChain`] when a face has three picked neighbours or is not joined
/// to the chain; [`GrooveRefusal::ArrowFaceNotAtToe`] when the outer face meets no picked face,
/// more than one, or one that is not an end of the chain. The errors of the ops calls otherwise.
pub fn chain(pick: &GroovePick<'_>, part: JointPart) -> Result<Chain> {
    let mut faces = pick.faces.to_vec();
    faces.sort_unstable();
    faces.dedup();
    if faces.is_empty() {
        return Err(GrooveRefusal::Invalid("no groove face picked"));
    }
    if faces.iter().chain([&pick.outer]).any(|&f| f >= pick.order.faces.len()) {
        return Err(GrooveRefusal::Invalid("a groove or outer face the part does not have"));
    }
    if faces.contains(&pick.outer) {
        return Err(GrooveRefusal::TJoint { part });
    }
    // each face's edges, the outer face's last
    let edges = faces.iter().chain([&pick.outer]).map(|&f| pick.part.rim_edges(pick.order, &[f])).collect::<solidcore_ops::Result<Vec<_>>>()?;
    let outer_edges = edges.last().map_or(&[][..], Vec::as_slice);
    let shared = |a: &[usize], b: &[usize]| a.iter().copied().filter(|e| b.contains(e)).collect::<Vec<usize>>();
    // each picked face's picked neighbours, by position in `faces`
    let mut next: Vec<Vec<usize>> = Vec::with_capacity(faces.len());
    for (i, f) in faces.iter().enumerate() {
        let near: Vec<usize> = (0..faces.len()).filter(|&j| j != i && !shared(&edges[i], &edges[j]).is_empty()).collect();
        if near.len() > 2 {
            return Err(GrooveRefusal::PicksNotOneChain { part, face: *f });
        }
        next.push(near);
    }
    // the toe: the one picked face beside the outer face, at an end of the chain (a closed ring
    // of faces has no end, so it lands here too)
    let beside: Vec<usize> = (0..faces.len()).filter(|&i| !shared(&edges[i], outer_edges).is_empty()).collect();
    let start = match beside.as_slice() {
        [i] if next[*i].len() <= 1 => *i,
        _ => return Err(GrooveRefusal::ArrowFaceNotAtToe { part, face: pick.outer }),
    };
    let mut order = vec![start];
    let mut at = start;
    while let Some(&j) = next[at].iter().find(|&&j| !order.contains(&j)) {
        order.push(j);
        at = j;
    }
    if let Some(off) = (0..faces.len()).find(|i| !order.contains(i)) {
        return Err(GrooveRefusal::PicksNotOneChain { part, face: faces[off] });
    }
    Ok(Chain { toe: shared(&edges[start], outer_edges), faces: order.into_iter().map(|i| faces[i]).collect() })
}

/// The joint's frame, from both picks and their chains: a plane whose normal is the joint's
/// direction `d`, along which every groove and outer face of both parts runs straight.
///
/// `d` is square to the two plane normals, among those faces, furthest from parallel (the
/// smallest `|n · m|`; the faces taken left part first, each part's groove faces from the toe,
/// then its outer face, and a tie keeps the first pair, so the order of the picks changes
/// nothing). Every one of those faces is then proved straight along `d`: it strays across it by
/// no more than the zone tolerance ([`Body::face_straying`]).
///
/// The frame (brett19 on #18500): `y` is the left outer face's outward normal with its part along
/// `d` taken out, `x = y × d`, and `d` takes the sign that makes `x` point from the left outer
/// face's centroid towards its toe, out of the left part; so swapping the parts reverses `d`. The
/// origin is the middle of the left toe's first edge. Every vector is an ops call
/// ([`affine::orthonormal_axes`], [`affine::rotation`], [`affine::distance_along`]).
///
/// # Errors
/// [`GrooveRefusal::OuterFaceNotPlane`]; [`GrooveRefusal::ArrowFaceNotAtToe`] when a part's
/// first groove face is a plane parallel to its outer face (they meet at no angle) or a chain has
/// no toe edge; [`GrooveRefusal::GrooveFaceNotStraight`]; [`GrooveRefusal::Invalid`] for a chain
/// with no face. The errors of the ops calls otherwise.
pub fn joint_frame(picks: [&GroovePick<'_>; 2], chains: [&Chain; 2]) -> Result<Plane> {
    Ok(frame_along(picks, chains)?.0)
}

/// [`joint_frame`] and the joint's direction, the plane's normal (`x × y`).
fn frame_along(picks: [&GroovePick<'_>; 2], chains: [&Chain; 2]) -> Result<(Plane, [f64; 3])> {
    let zone = sweep::zone_tolerance()?;
    let parts = [JointPart::Left, JointPart::Right];
    let mut outer = [[0.0; 3]; 2];
    // every face of the joint: left part first, groove faces from the toe, then the outer face
    let mut faces: Vec<(JointPart, &GroovePick<'_>, usize)> = Vec::new();
    for (((pick, chain), part), out) in picks.into_iter().zip(chains).zip(parts).zip(&mut outer) {
        let Some(n) = pick.part.face_plane_normal(pick.order, pick.outer)? else {
            return Err(GrooveRefusal::OuterFaceNotPlane { part, face: pick.outer });
        };
        *out = n.value;
        let no_toe = GrooveRefusal::ArrowFaceNotAtToe { part, face: pick.outer };
        let Some(&first) = chain.faces.first() else {
            return Err(GrooveRefusal::Invalid("a chain with no groove face"));
        };
        // a first groove face parallel to the outer face meets it at no angle: there is no toe
        let inside = pick.part.face_centroid(pick.order, pick.outer, 1e-8)?.centroid;
        if chain.toe.is_empty() || pick.part.face_square_level(pick.order, first, n.value, inside)?.is_some() {
            return Err(no_toe);
        }
        faces.extend(chain.faces.iter().chain([&pick.outer]).map(|&f| (part, pick, f)));
    }

    // the direction: square to the two plane normals furthest from parallel
    let mut normals = Vec::new();
    for &(_, pick, f) in &faces {
        if let Some(n) = pick.part.face_plane_normal(pick.order, f)? {
            normals.push(n.value);
        }
    }
    let mut best: Option<(f64, [f64; 3], [f64; 3])> = None;
    for (i, &n) in normals.iter().enumerate() {
        for &m in normals.iter().skip(i + 1) {
            let c = affine::distance_along([0.0; 3], n, m).abs();
            if best.is_none_or(|(b, ..)| c < b) {
                best = Some((c, n, m));
            }
        }
    }
    // the two outer faces are planes, so a pair is always there
    let Some((_, n, m)) = best else {
        return Err(GrooveRefusal::Invalid("fewer than two plane faces"));
    };
    // d = n × m', m' the unit part of m square to n: m' turned a quarter about n
    let [n, m] = affine::orthonormal_axes(n, m)?;
    let d = affine::apply_point(&affine::rotation([0.0; 3], n, FRAC_PI_2)?, m);
    for &(part, pick, f) in &faces {
        if !pick.part.face_straying(pick.order, f, d)?.is_some_and(|e| e <= zone) {
            return Err(GrooveRefusal::GrooveFaceNotStraight { part, face: f });
        }
    }

    // the frame, from the left part: y the outer normal square to d, x = y × d (d turned a quarter
    // about y), its sign towards the toe
    let [left, _] = picks;
    let [_, y] = affine::orthonormal_axes(d, outer[0])?;
    let mut x = affine::apply_point(&affine::rotation([0.0; 3], y, FRAC_PI_2)?, d);
    let Some(&toe_edge) = chains[0].toe.first() else {
        return Err(GrooveRefusal::ArrowFaceNotAtToe { part: JointPart::Left, face: left.outer });
    };
    let toe = left.part.edge_middle(left.order, toe_edge)?;
    let inside = left.part.face_centroid(left.order, left.outer, 1e-8)?.centroid;
    let mut d = d;
    if affine::distance_along(inside, x, toe) < 0.0 {
        x = x.map(|v| -v);
        d = d.map(|v| -v);
    }
    let v = |[a, b, c]: [f64; 3]| Vec3::new(a, b, c);
    Ok((Plane { origin: v(toe), x: v(x), y: v(y) }, d))
}

/// A groove weld: its body and what built it.
#[derive(Debug)]
pub struct GrooveWeld {
    /// The weld's body: the fill swept along the joint over the parts' length.
    pub body: Body,
    /// The joint's frame ([`joint_frame`]): the section lies in it, the body runs along its normal.
    pub frame: Plane,
    /// The section's fill: its region in the frame, its depths, what a caller may report.
    pub fill: ChannelFill,
}

/// The weld that fills the groove between the two picks, as `options` asks: each part's chain
/// ([`chain`]), the joint's frame ([`joint_frame`]), each part's wall read in that frame from its
/// groove faces' edges at the ends of the joint ([`planar_edges`]), the channel between the walls
/// filled ([`fill_channel`]), the fill swept along the joint over the parts' length
/// ([`extrude_planar_regions`]). The depth is measured from the left outer face's level where
/// both outer faces are square to the frame's `y` (a butt joint), else from the chord between
/// the toes. Nothing writes to a part.
///
/// # Errors
/// The refusals of [`chain`] and [`joint_frame`]; [`GrooveRefusal::SectionChanges`];
/// [`GrooveRefusal::PartsNotAlongside`]; [`GrooveRefusal::Fill`] with the size that builds. The
/// errors of the ops calls otherwise, a fill's input the kernel refuses among them.
pub fn groove_weld(picks: [&GroovePick<'_>; 2], options: &FillOptions) -> Result<GrooveWeld> {
    let parts = [JointPart::Left, JointPart::Right];
    let chains = [chain(picks[0], parts[0])?, chain(picks[1], parts[1])?];
    let (frame, d) = frame_along(picks, [&chains[0], &chains[1]])?;
    let zone = sweep::zone_tolerance()?;
    let mut walls: [Vec<Seg2>; 2] = [Vec::new(), Vec::new()];
    for (((pick, chain), part), wall) in picks.into_iter().zip(&chains).zip(parts).zip(&mut walls) {
        *wall = wall_of(pick, chain, part, &frame, zone)?;
    }

    // both toes the same length along the joint, side by side
    let mut toes = [(0.0, [0.0; 3]); 2];
    for ((pick, chain), toe) in picks.into_iter().zip(&chains).zip(&mut toes) {
        let Some(&edge) = chain.toe.first() else {
            return Err(GrooveRefusal::Invalid("a chain with no toe edge"));
        };
        *toe = (pick.part.edge_length(pick.order, edge, 1e-9)?, pick.part.edge_middle(pick.order, edge)?);
    }
    let [(length, _), (other, middle)] = toes;
    let o = [frame.origin.x, frame.origin.y, frame.origin.z];
    if (length - other).abs() > zone || affine::distance_along(o, d, middle).abs() > zone {
        return Err(GrooveRefusal::PartsNotAlongside);
    }

    // a butt joint measures the depth from the plates' outer faces, square to the frame's y
    let y = [frame.y.x, frame.y.y, frame.y.z];
    let right = picks[1];
    let inside = right.part.face_centroid(right.order, right.outer, 1e-8)?.centroid;
    let butt = right.part.face_square_level(right.order, right.outer, y, inside)?.is_some();
    let depth_normal = butt.then_some(Vec2::new(0.0, 1.0));
    let [left_wall, right_wall] = &walls;
    let fill = match fill_channel(left_wall, right_wall, depth_normal, options)? {
        Ok(fill) => fill,
        Err(refusal) => return Err(GrooveRefusal::Fill { why: refusal.why, fit: refusal.fit }),
    };
    let body = extrude_planar_regions(&frame, std::slice::from_ref(&fill.region), -length / 2.0, length / 2.0)?;
    Ok(GrooveWeld { body, frame, fill })
}

/// A part's wall in the joint's frame, from its toe away from the face: each groove face's
/// section, read from its edges at the ends of the joint, which must agree, each starting where
/// the one before ends.
fn wall_of(pick: &GroovePick<'_>, chain: &Chain, part: JointPart, frame: &Plane, zone: f64) -> Result<Vec<Seg2>> {
    let section = |face: usize| -> Result<Seg2> {
        let edges = pick.part.rim_edges(pick.order, &[face])?;
        let read = planar_edges(pick.part, pick.order, frame, &edges)?;
        let tol = read.tolerance.unwrap_or(0.0).max(zone);
        let ends: Vec<Seg2> = read.pieces.into_iter().flatten().filter_map(|p| one_seg(&p)).collect();
        match ends.as_slice() {
            [a, b] if same(a, b, tol) => Ok(*a),
            _ => Err(GrooveRefusal::SectionChanges { part, face }),
        }
    };
    let tol = zone * 8.0;
    let outer = section(pick.outer)?;
    let mut wall: Vec<Seg2> = Vec::new();
    for &face in &chain.faces {
        let seg = section(face)?;
        let (a, b) = ends(&seg);
        // the first face starts at the toe, where it meets the outer face
        let at = wall.last().map_or_else(|| ends_meeting(&outer, a, b, tol), |last| Some(ends(last).1));
        let seg = match at {
            Some(p) if near(a, p, tol) => seg,
            Some(p) if near(b, p, tol) => reversed(&seg),
            _ => return Err(GrooveRefusal::PicksNotOneChain { part, face }),
        };
        wall.push(seg);
    }
    Ok(wall)
}

/// The piece of one edge's read when it is a single line or arc.
fn one_seg(pieces: &[PlanarPiece]) -> Option<Seg2> {
    match pieces {
        [PlanarPiece::Seg(s @ (Seg2::Line { .. } | Seg2::Arc { .. }))] => Some(*s),
        _ => None,
    }
}

/// A line's or an arc's two ends. shortcut: an arc's ends by its angles in plain trigonometry,
/// as the app's sheet-metal code reads them; an ops call for it if reviewers ask.
fn ends(s: &Seg2) -> (Vec2, Vec2) {
    match *s {
        Seg2::Arc { center, radius, start, sweep } => {
            let at = |t: f64| Vec2::new(center.x + radius * t.cos(), center.y + radius * t.sin());
            (at(start), at(start + sweep))
        }
        Seg2::Line { a, b } => (a, b),
        _ => (Vec2::new(f64::NAN, f64::NAN), Vec2::new(f64::NAN, f64::NAN)),
    }
}

fn near(p: Vec2, q: Vec2, tol: f64) -> bool {
    (p.x - q.x).abs() <= tol && (p.y - q.y).abs() <= tol
}

/// Which of `a` and `b` is an end of `s`.
fn ends_meeting(s: &Seg2, a: Vec2, b: Vec2, tol: f64) -> Option<Vec2> {
    let (p, q) = ends(s);
    [a, b].into_iter().find(|&e| near(e, p, tol) || near(e, q, tol))
}

/// Whether two reads of a face's ends are the same line or arc, either way round.
fn same(a: &Seg2, b: &Seg2, tol: f64) -> bool {
    let ((a0, a1), (b0, b1)) = (ends(a), ends(b));
    let ends_agree = (near(a0, b0, tol) && near(a1, b1, tol)) || (near(a0, b1, tol) && near(a1, b0, tol));
    let kind_agrees = match (a, b) {
        (Seg2::Line { .. }, Seg2::Line { .. }) => true,
        (Seg2::Arc { center: c, radius: r, .. }, Seg2::Arc { center: k, radius: q, .. }) => near(*c, *k, tol) && (r - q).abs() <= tol,
        _ => false,
    };
    ends_agree && kind_agrees
}

/// The same line or arc the other way round.
fn reversed(s: &Seg2) -> Seg2 {
    match *s {
        Seg2::Line { a, b } => Seg2::Line { a: b, b: a },
        Seg2::Arc { center, radius, start, sweep } => Seg2::Arc { center, radius, start: start + sweep, sweep: -sweep },
        other => other,
    }
}

#[cfg(test)]
#[path = "weld_groove_tests.rs"]
#[allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]
mod tests;
