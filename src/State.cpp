#include "generalized_information_filter/State.hpp"

namespace GIF {

StateBase::StateBase(const std::shared_ptr<const StateDefinition>& def)
    : def_(def) {
}

StateBase::~StateBase() {}

bool StateBase::matchesDef(
    const std::shared_ptr<const StateDefinition>& def) const {
  if (getNumElement() != def->getNumElement()) {
    return false;
  }
  for (int i = 0; i < getNumElement(); i++) {
    if (!def->getElementDefinition(i)->isOfDef(getElement(i))) {
      return false;
    }
  }
  return true;
}

StateBase& StateBase::operator=(const StateBase& other) {
  for (int i = 0; i < getNumElement(); i++) {
    *getElement(i) = *other.getElement(i);
  }
  return *this;
}

void StateBase::print() const {
  for (int i = 0; i < getNumElement(); i++) {
    std::cout << getDef()->getName(i) << ": ";
    getElement(i)->print();
  }
}

void StateBase::setIdentity() {
  for (int i = 0; i < getNumElement(); i++) {
    getElement(i)->setIdentity();
  }
}

void StateBase::setRandom(int& s) {
  for (int i = 0; i < getNumElement(); i++) {
    getElement(i)->setRandom(s);
  }
}

void StateBase::boxplus(const Eigen::Ref<const Eigen::VectorXd>& vec,
                        const std::shared_ptr<StateBase>& out) const {
  for (int i = 0; i < getNumElement(); i++) {
    getElement(i)->boxplus(
        vec.block(getStart(i), 0, getElement(i)->getDim(), 1),
        out->getElement(i));
  }
}

void StateBase::boxminus(const std::shared_ptr<const StateBase>& ref,
                         Eigen::Ref<Eigen::VectorXd> vec) const {
  for (int i = 0; i < getNumElement(); i++) {
    getElement(i)->boxminus(
        ref->getElement(i),
        vec.block(getStart(i), 0, getElement(i)->getDim(), 1));
  }
}

std::shared_ptr<const StateDefinition> StateBase::getDef() const {
  return def_;
}

State::State(const std::shared_ptr<const StateDefinition>& def)
    : StateBase(def) {
  for (int i = 0; i < def_->getNumElement(); i++) {
    elements_.push_back(def_->getElementDefinition(i)->newElement());
  }
}
;

State::~State() {
}
;
State& State::operator=(const StateBase& other) {
  dynamic_cast<const StateBase&>(*this) = other;
  return *this;
}

int State::getNumElement() const {
  return elements_.size();
}

StateWrapper::StateWrapper(const std::shared_ptr<const StateDefinition>& def,
                           const std::shared_ptr<const StateDefinition>& in)
    : StateBase(def),
      in_(in) {
  computeMap();
}
;

StateWrapper::~StateWrapper() {
}
;
StateWrapper& StateWrapper::operator=(const StateBase& other) {
  dynamic_cast<const StateBase&>(*this) = other;
  return *this;
}

int StateWrapper::getNumElement() const {
  return indexMap_.size();
}

void StateWrapper::computeMap() {
  indexMap_.resize(def_->getNumElement());
  for (int i = 0; i < def_->getNumElement(); i++) {
    indexMap_[i] = in_->findName(def_->getName(i));
    assert(indexMap_.at(i) != -1);
  }
}

void StateWrapper::setState(const std::shared_ptr<StateBase>& state) {
  state->matchesDef(in_);
  state_ = state;
  constState_ = state;
}

void StateWrapper::setState(
    const std::shared_ptr<const StateBase>& state) const {
  state->matchesDef(in_);
  constState_ = state;
}

void StateWrapper::wrapJacobian(Eigen::Ref<MXD> out,
                                const Eigen::Ref<const MXD>& in,
                                int rowOffset) const {
  const int rows = in.rows();
  for (int i = 0; i < getNumElement(); ++i) {
    const int cols = getElement(i)->getDim();
    out.block(rowOffset, constState_->getStart(indexMap_.at(i)), rows, cols) =
        in.block(0, getStart(i), rows, cols);
  }
}

}
