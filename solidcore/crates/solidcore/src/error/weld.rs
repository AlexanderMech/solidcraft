//! Why a weld is refused, in solidcore's terms: one enum for every weld type (`crate::weld`).
//!
//! Not upstream: SolidCraft's own (#15132, #18500). A weld is display geometry beside the
//! parts it joins; these say why the weld the user defined cannot be built as defined. Where
//! a smaller weld builds, the failure's [`Fit`](super::Fit) names it. Bad input (walls that are
//! not lines and arcs, sizes that are not positive) stays
//! [`SolidError::InvalidInput`](super::SolidError::InvalidInput).
//!
//! | Rust                                      | Status difftest compares |
//! |-------------------------------------------|--------------------------|
//! | [`WeldFailure::ConcaveFaceCrossesWall`]   | `SM_ERR` (no C++)        |
//! | [`WeldFailure::ConcaveFaceTooDeep`]       | `SM_ERR` (no C++)        |
//! | [`WeldFailure::ConvexFaceTooHigh`]        | `SM_ERR` (no C++)        |
//! | [`WeldFailure::BlendedToesOffPlane`]      | `SM_ERR` (no C++)        |

use super::{ErrorKind, SM_ERR};

/// Which side of a weld something is about.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum WeldSide {
    /// The face (the side the weld is made from; the side its walls start from).
    Face,
    /// The second face of a weld from both sides (a weld symbol above and below the
    /// reference line).
    Second,
}

/// Why a weld cannot be built as defined.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum WeldFailure {
    /// A concave face `depth` deep reaches the section's other sides (a groove wall, the other
    /// face or the bottom). The failure's fit is the depth where the face first touches one,
    /// in closed form; every depth below it less the kernel's zone tolerance builds.
    ConcaveFaceCrossesWall {
        /// Which face.
        side: WeldSide,
        /// The depth asked for.
        depth: f64,
    },
    /// A concave face `depth` deep dips past half a circle, touching nothing. The failure's
    /// fit is the half circle: half the toe distance, which builds.
    ConcaveFaceTooDeep {
        /// Which face.
        side: WeldSide,
        /// The depth asked for.
        depth: f64,
    },
    /// A convex face `height` high is past half a circle (its main arc's half sweep past a
    /// quarter turn). The failure's fit is the highest that builds: half the toe distance,
    /// plus the blend's overlap for blended toes.
    ConvexFaceTooHigh {
        /// Which face.
        side: WeldSide,
        /// The height asked for.
        height: f64,
    },
    /// Blended toes where the plates' outer faces on that side do not run along the chord
    /// between the toes (offset plates, plates of unequal thickness on that side, or no
    /// common plate normal, as in a T- or corner joint): the run-out would not be tangent to
    /// either plate.
    BlendedToesOffPlane {
        /// Which face.
        side: WeldSide,
    },
}

impl WeldFailure {
    /// The status the C++ reports for it (a comparison detail): `SM_ERR`; the C++ has no
    /// welds.
    pub(crate) const fn comparison_status(&self) -> i32 {
        match self {
            Self::ConcaveFaceCrossesWall { .. }
            | Self::ConcaveFaceTooDeep { .. }
            | Self::ConvexFaceTooHigh { .. }
            | Self::BlendedToesOffPlane { .. } => SM_ERR,
        }
    }

    /// Its stable [`ErrorKind`]: a size or a combination the weld does not support.
    #[must_use]
    pub const fn kind(&self) -> ErrorKind {
        match self {
            Self::ConcaveFaceCrossesWall { .. }
            | Self::ConcaveFaceTooDeep { .. }
            | Self::ConvexFaceTooHigh { .. }
            | Self::BlendedToesOffPlane { .. } => ErrorKind::Unsupported,
        }
    }
}

impl core::fmt::Display for WeldFailure {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        let name = |side: WeldSide| match side {
            WeldSide::Face => "face",
            WeldSide::Second => "second face",
        };
        match *self {
            Self::ConcaveFaceCrossesWall { side, depth } => write!(
                f,
                "a concave {} {depth} deep reaches the section's other sides",
                name(side)
            ),
            Self::ConcaveFaceTooDeep { side, depth } => write!(
                f,
                "a concave {} {depth} deep is past half a circle",
                name(side)
            ),
            Self::ConvexFaceTooHigh { side, height } => write!(
                f,
                "a convex {} {height} high is past half a circle",
                name(side)
            ),
            Self::BlendedToesOffPlane { side } => write!(
                f,
                "blended toes on a {} whose plates do not run along its chord",
                name(side)
            ),
        }
    }
}
