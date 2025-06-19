use crate::asthenosphere::AsthenosphereCell;
use std::cell::RefCell;
use std::rc::{Rc, Weak};
use h3o::CellIndex;

#[derive(Debug)]
pub struct AsthenosphereCellLinked {
    pub cell: AsthenosphereCell,
    pub next: Option<Rc<RefCell<AsthenosphereCellLinked>>>,
    pub prev: Option<Weak<RefCell<AsthenosphereCellLinked>>>,
}

impl AsthenosphereCellLinked {
    pub fn from_cell(cell: &AsthenosphereCell) -> AsthenosphereCellLinked {
        AsthenosphereCellLinked {
            next: None,
            prev: None,
            cell: cell.clone()
        }
    }

    pub fn new(cell: AsthenosphereCell) -> Self {
        AsthenosphereCellLinked {
            cell,
            next: None,
            prev: None,
        }
    }

    pub fn add(self_rc: &Rc<RefCell<AsthenosphereCellLinked>>) -> Rc<RefCell<Self>> {
        if let Some(next_rc) = &self_rc.borrow().next {
            Rc::clone(next_rc);
        }

        let new_cell = AsthenosphereCell {
            step: self_rc.borrow().cell.step + 1,
            ..self_rc.borrow().cell.clone()
        };

        let new_linked = Rc::new(RefCell::new(AsthenosphereCellLinked::new(new_cell)));

        self_rc.borrow_mut().next = Some(Rc::clone(&new_linked));
        new_linked.borrow_mut().prev = Some(Rc::downgrade(self_rc));

        new_linked
    }

    pub fn unlink_next(&mut self) {
        if let Some(next_rc) = self.next.take() {
            next_rc.borrow_mut().prev = None;
        }
    }

    pub fn unlink_all_prev(&mut self) {
        if let Some(prev_weak) = self.prev.take() {
            if let Some(prev_rc) = prev_weak.upgrade() {
                prev_rc.borrow_mut().unlink_all_prev();
                prev_rc.borrow_mut().next = None;
            }
        }
    }

    pub fn unlink_all_next(&mut self) {
        if let Some(next_rc) = self.next.take() {
            next_rc.borrow_mut().unlink_all_next();
            next_rc.borrow_mut().prev = None;
        }
    }
}
